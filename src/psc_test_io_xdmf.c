
#define BOUNDS_CHECK
#define BOUNDSCHECK

#include "psc_test_io_xdmf.h"

#include <mrc_domain.h>

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// ----------------------------------------------------------------------
// xdmf_collective_setup

void
xdmf_collective_setup(struct xdmf *xdmf)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; MPI_Comm_rank(comm, &rank);
  int size; MPI_Comm_size(comm, &size);

  if (rank < xdmf->nr_writers) {
    xdmf->is_writer = 1;
  }
  MPI_Comm_split(comm, xdmf->is_writer, rank, &xdmf->comm_writers);
}

// ----------------------------------------------------------------------
// xdmf_collective_destroy

void
xdmf_collective_destroy(struct xdmf *xdmf)
{
  if (xdmf->comm_writers) {
    MPI_Comm_free(&xdmf->comm_writers);
  }
}

// ======================================================================

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static bool
find_intersection(int *ilo, int *ihi, const int *ib1, const int *im1,
		  const int *ib2, const int *im2)
{
  bool has_intersection = true;
  for (int d = 0; d < 3; d++) {
    ilo[d] = MAX(ib1[d], ib2[d]);
    ihi[d] = MIN(ib1[d] + im1[d], ib2[d] + im2[d]);
    if (ihi[d] - ilo[d] <= 0) {
      has_intersection = false;
    }
  }
  return has_intersection;
}

// ----------------------------------------------------------------------
// collective helper context

struct collective_m3_recv_patch {
  int ilo[3]; // intersection low
  int ihi[3]; // intersection high
};

struct collective_m3_peer {
  int rank;
  struct collective_m3_recv_patch *begin;
  struct collective_m3_recv_patch *end;
  float *recv_buf;
  size_t buf_size;
};

struct collective_m3_ctx {
  int gdims[3];
  int slab_dims[3], slab_off[3];
  int nr_patches, nr_global_patches;
  int slow_dim;
  int slow_indices_per_writer;
  int slow_indices_rmndr;

  struct collective_m3_recv_patch *recv_patches;
  int n_recv_patches;

  MPI_Request *recv_reqs;
  struct collective_m3_peer *peers;
  int n_peers;

  float **send_bufs; // one for each writer
  MPI_Request *send_reqs;
  int nr_sends;

  MPI_Comm comm;
  int rank;
  int size;
};

static void
collective_m3_init(struct xdmf *xdmf, struct collective_m3_ctx *ctx,
		   struct mrc_domain *domain)
{
  mrc_domain_get_global_dims(domain, ctx->gdims);
  mrc_domain_get_patches(domain, &ctx->nr_patches);
  mrc_domain_get_nr_global_patches(domain, &ctx->nr_global_patches);
  for (int d = 0; d < 3; d++) {
    ctx->slab_dims[d] = ctx->gdims[d];
    ctx->slab_off[d] = 0;
  }
  ctx->slow_dim = 2;
  while (ctx->gdims[ctx->slow_dim] == 1) {
    ctx->slow_dim--;
  }
  assert(ctx->slow_dim >= 0);
  int total_slow_indices = ctx->slab_dims[ctx->slow_dim];
  ctx->slow_indices_per_writer = total_slow_indices / xdmf->nr_writers;
  ctx->slow_indices_rmndr = total_slow_indices % xdmf->nr_writers;

  ctx->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(ctx->comm, &ctx->rank);
  MPI_Comm_size(ctx->comm, &ctx->size);
}

static void
get_writer_off_dims(struct collective_m3_ctx *ctx, int writer,
		    int *writer_off, int *writer_dims)
{
  for (int d = 0; d < 3; d++) {
    writer_dims[d] = ctx->slab_dims[d];
    writer_off[d] = ctx->slab_off[d];
  }
  writer_dims[ctx->slow_dim] = ctx->slow_indices_per_writer + (writer < ctx->slow_indices_rmndr);
  if (writer < ctx->slow_indices_rmndr) {
    writer_off[ctx->slow_dim] += (ctx->slow_indices_per_writer + 1) * writer;
  } else {
    writer_off[ctx->slow_dim] += ctx->slow_indices_rmndr +
      ctx->slow_indices_per_writer * writer;
  }
}

// ----------------------------------------------------------------------
// collective_send_fld_begin

static void
collective_send_fld_begin(struct xdmf *xdmf, struct collective_m3_ctx *ctx,
			  struct mrc_domain *domain, int m)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);

  size_t *buf_sizes = calloc(xdmf->nr_writers, sizeof(*buf_sizes));
  ctx->send_bufs = calloc(xdmf->nr_writers, sizeof(*ctx->send_bufs));
  ctx->send_reqs = calloc(xdmf->nr_writers, sizeof(*ctx->send_reqs));

  for (int writer = 0; writer < xdmf->nr_writers; writer++) {
    // don't send to self
    if (writer == ctx->rank) {
      ctx->send_reqs[writer] = MPI_REQUEST_NULL;
      continue;
    }
    int writer_off[3], writer_dims[3];
    get_writer_off_dims(ctx, writer, writer_off, writer_dims);

    // find buf_size per writer
    for (int p = 0; p < nr_patches; p++) {
      int *ldims = patches[p].ldims;
      int ilo[3], ihi[3];
      bool has_intersection =
	find_intersection(ilo, ihi, patches[p].off, ldims,
			  writer_off, writer_dims);
      if (!has_intersection)
	continue;

      buf_sizes[writer] += ldims[0] * ldims[1] * ldims[2];
    }

    // allocate buf per writer
    //mprintf("to writer %d buf_size %d\n", writer, buf_sizes[writer]);
    if (buf_sizes[writer] == 0) {
      ctx->send_reqs[writer] = MPI_REQUEST_NULL;
      continue;
    }

    ctx->send_bufs[writer] = malloc(buf_sizes[writer] * sizeof(*ctx->send_bufs[writer]));
    assert(ctx->send_bufs[writer]);
    
    MPI_Isend(ctx->send_bufs[writer], buf_sizes[writer], MPI_FLOAT,
	      writer, 0x1000, ctx->comm,
	      &ctx->send_reqs[writer]);
  }
  free(buf_sizes);
}

// ----------------------------------------------------------------------
// collective_send_fld_end

static void
collective_send_fld_end(struct xdmf *xdmf, struct collective_m3_ctx *ctx)
{
  MPI_Waitall(xdmf->nr_writers, ctx->send_reqs, MPI_STATUSES_IGNORE);

  for (int writer = 0; writer < xdmf->nr_writers; writer++) {
    free(ctx->send_bufs[writer]);
  }
  free(ctx->send_bufs);
  free(ctx->send_reqs);
}
    
// ----------------------------------------------------------------------
// writer_comm_init

static void
writer_comm_init(struct collective_m3_ctx *ctx, int *writer_offs, int *writer_dims,
		 struct mrc_domain *domain, int size_of_type)
{
  // find out who's sending, OPT: this way is not very scalable
  // could also be optimized by just looking at slow_dim
  // FIXME, figure out pattern and cache, at least across components

  int nr_global_patches;
  mrc_domain_get_nr_global_patches(domain, &nr_global_patches);

  ctx->n_recv_patches = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);

    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     writer_offs, writer_dims);
    if (!has_intersection) {
      continue;
    }

    ctx->n_recv_patches++;
  }
  mprintf("n_recv_patches %d\n", ctx->n_recv_patches);

  ctx->recv_patches = calloc(ctx->n_recv_patches, sizeof(*ctx->recv_patches));

  struct collective_m3_recv_patch **recv_patches_by_rank = calloc(ctx->size + 1, sizeof(*recv_patches_by_rank));

  int cur_rank = -1;
  ctx->n_recv_patches = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);

    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     writer_offs, writer_dims);
    if (!has_intersection) {
      continue;
    }

    assert(info.rank >= cur_rank);
    while (cur_rank < info.rank) {
      cur_rank++;
      recv_patches_by_rank[cur_rank] = &ctx->recv_patches[ctx->n_recv_patches];
      //mprintf("rank %d patches start at %d\n", cur_rank, ctx->n_recv_patches);
    }

    struct collective_m3_recv_patch *recv_patch = &ctx->recv_patches[ctx->n_recv_patches];
    for (int d = 0; d < 3; d++) {
      recv_patch->ilo[d] = ilo[d];
      recv_patch->ihi[d] = ihi[d];
    }
    ctx->n_recv_patches++;
  }

  while (cur_rank < ctx->size) {
    cur_rank++;
    recv_patches_by_rank[cur_rank] = &ctx->recv_patches[ctx->n_recv_patches];
    //mprintf("rank %d patches start at %d\n", cur_rank, ctx->n_recv_patches);
  }

  ctx->n_peers = 0;
  for (int rank = 0; rank < ctx->size; rank++) {
    struct collective_m3_recv_patch *begin = recv_patches_by_rank[rank];
    struct collective_m3_recv_patch *end   = recv_patches_by_rank[rank+1];

    if (begin == end) {
      continue;
    }

    //mprintf("peer rank %d # = %ld\n", rank, end - begin);
    ctx->n_peers++;
  }
  mprintf("n_peers %d\n", ctx->n_peers);

  ctx->peers = calloc(ctx->n_peers, sizeof(*ctx->peers));
  ctx->n_peers = 0;
  for (int rank = 0; rank < ctx->size; rank++) {
    struct collective_m3_recv_patch *begin = recv_patches_by_rank[rank];
    struct collective_m3_recv_patch *end   = recv_patches_by_rank[rank+1];

    if (begin == end) {
      continue;
    }

    struct collective_m3_peer *peer = &ctx->peers[ctx->n_peers];
    peer->rank = rank;
    peer->begin = begin;
    peer->end = end;

    // for remote patches, allocate buffer
    if (peer->rank != ctx->rank) {
      peer->buf_size = 0;
      for (struct collective_m3_recv_patch *recv_patch = peer->begin; recv_patch < peer->end; recv_patch++) {
	int *ilo = recv_patch->ilo, *ihi = recv_patch->ihi;
	peer->buf_size += (size_t) (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);
      }
      
      // alloc aggregate recv buffers
      peer->recv_buf = malloc(peer->buf_size * sizeof(*peer->recv_buf));
    }

    ctx->n_peers++;
  }
  free(recv_patches_by_rank);

  ctx->recv_reqs = calloc(ctx->n_peers, sizeof(*ctx->recv_reqs));
  mprintf("recv_reqs %p\n", ctx->recv_reqs);
}

// ----------------------------------------------------------------------
// writer_comm_begin

static void
writer_comm_begin(struct collective_m3_ctx *ctx)
{
  for (struct collective_m3_peer *peer = ctx->peers; peer < ctx->peers + ctx->n_peers; peer++) {
    // skip local patches
    if (peer->rank == ctx->rank) {
      ctx->recv_reqs[peer - ctx->peers] = MPI_REQUEST_NULL;
      continue;
    }

    // recv aggregate buffers
    MPI_Irecv(peer->recv_buf, peer->buf_size, MPI_FLOAT, peer->rank, 0x1000, ctx->comm,
	      &ctx->recv_reqs[peer - ctx->peers]);
  }
}

// ----------------------------------------------------------------------
// writer_comm_end

static void
writer_comm_end(struct collective_m3_ctx *ctx)
{
  MHERE;
  MPI_Waitall(ctx->n_peers, ctx->recv_reqs, MPI_STATUSES_IGNORE);
  MHERE;
  
  for (struct collective_m3_peer *peer = ctx->peers; peer < ctx->peers + ctx->n_peers; peer++) {
    // skip local patches
    if (peer->rank == ctx->rank) {
      continue;
    }

    // copy
  }
}

// ----------------------------------------------------------------------
// writer_comm_destroy

static void
writer_comm_destroy(struct collective_m3_ctx *ctx)
{
  free(ctx->recv_reqs);
  
  for (struct collective_m3_peer *peer = ctx->peers; peer < ctx->peers + ctx->n_peers; peer++) {
    free(peer->recv_buf);
  }
  free(ctx->peers);

  free(ctx->recv_patches);
}

// ----------------------------------------------------------------------
// xdmf_collective_write_m3

void
xdmf_collective_write_m3(struct xdmf* xdmf, const char *path, struct mrc_domain *domain)
{
  struct collective_m3_ctx ctx;
  collective_m3_init(xdmf, &ctx, domain);

  if (xdmf->is_writer) {
    int writer;
    MPI_Comm_rank(xdmf->comm_writers, &writer);
    int writer_dims[3], writer_off[3];
    get_writer_off_dims(&ctx, writer, writer_off, writer_dims);
    mprintf("writer_off %d %d %d dims %d %d %d\n",
    	    writer_off[0], writer_off[1], writer_off[2],
    	    writer_dims[0], writer_dims[1], writer_dims[2]);

    writer_comm_init(&ctx, writer_off, writer_dims, domain, sizeof(float));
    for (int m = 0; m < 2; m++) {
      writer_comm_begin(&ctx);
      collective_send_fld_begin(xdmf, &ctx, domain, 0);
      writer_comm_end(&ctx);
      collective_send_fld_end(xdmf, &ctx);
    }
    writer_comm_destroy(&ctx);
  } else {
    for (int m = 0; m < 2; m++) {
      collective_send_fld_begin(xdmf, &ctx, domain, 0);
      collective_send_fld_end(xdmf, &ctx);
    }
  }
}

