
#define DDCP_TYPE_COMMON     1
#define DDCP_TYPE_COMMON2    2
#define DDCP_TYPE_COMMON_OMP 3
#define DDCP_TYPE_CUDA       4

#include <mrc_ddc.h>

#include <string.h>

#define N_DIR (27)

// ----------------------------------------------------------------------
// ddcp_buf_t

typedef struct {
  int m_size;
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP
  struct psc_mparticles *m_mprts;
  int m_p;
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  struct cuda_bnd *m_bpatch;
#endif
} ddcp_buf_t;

static void
ddcp_buf_ctor(ddcp_buf_t *buf, struct psc_mparticles *mprts, int p)
{
  buf->m_size = 0;
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP
  buf->m_mprts = mprts;
  buf->m_p = p;
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  buf->m_bpatch = &cmprts->bnd.bpatch[p];
#endif
}

static void
ddcp_buf_dtor(ddcp_buf_t *buf)
{
}

static particle_t *
ddcp_buf_at(ddcp_buf_t *buf, int n)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP
  return mparticles_get_one(buf->m_mprts, buf->m_p, n);
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  return &buf->m_bpatch->prts[n];
#endif
}

static void
ddcp_buf_reserve(ddcp_buf_t *buf, int new_capacity)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP
  mparticles_patch_reserve(buf->m_mprts, buf->m_p, new_capacity);
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  if (buf->m_bpatch->capacity >= new_capacity) {
    return;
  }

  buf->m_bpatch->capacity = MAX(new_capacity, 2 * buf->m_bpatch->capacity);
  buf->m_bpatch->prts = realloc(buf->m_bpatch->prts,
				buf->m_bpatch->capacity * sizeof(*buf->m_bpatch->prts));
#endif
}

static unsigned int
ddcp_buf_capacity(ddcp_buf_t *buf)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP
  return mparticles_patch_capacity(buf->m_mprts, buf->m_p);
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  return buf->m_bpatch->capacity;
#endif
}

static unsigned int
ddcp_buf_size(ddcp_buf_t *buf)
{
  return buf->m_size;
}

static void
ddcp_buf_resize(ddcp_buf_t *buf, int new_size)
{
  assert(new_size <= ddcp_buf_capacity(buf));
  buf->m_size = new_size;
}

static void
ddcp_buf_push_back(ddcp_buf_t *buf, particle_t *p)
{
  // this assert should go away in the general case,
  // but we actually push_back into the same array we're reading from,
  // so we better don't go beyond current capacity
  assert(buf->m_size + 1 <= ddcp_buf_capacity(buf));
  ddcp_buf_reserve(buf, buf->m_size + 1);
  *ddcp_buf_at(buf, buf->m_size) = *p;
  buf->m_size++;
}

// ----------------------------------------------------------------------
// particle_buf_t

typedef struct {
  particle_t *m_data;
  unsigned int m_size;
  unsigned int m_capacity;
} particle_buf_t;

static void
particle_buf_ctor(particle_buf_t *buf)
{
  buf->m_capacity = 8;
  buf->m_data = malloc(buf->m_capacity * sizeof(*buf->m_data));
  buf->m_size = 0;
}

static void
particle_buf_dtor(particle_buf_t *buf)
{
  free(buf->m_data);
}

static unsigned int
particle_buf_size(particle_buf_t *buf)
{
  return buf->m_size;
}

static void
particle_buf_resize(particle_buf_t *buf, unsigned int size)
{
  assert(size <= buf->m_capacity);
  buf->m_size = size;
}

static void
particle_buf_reserve(particle_buf_t *buf, unsigned int new_capacity)
{
  if (new_capacity > buf->m_capacity) {
    // reallocate a larger buffer, at least doubling buffer size each time
    buf->m_capacity = MAX(new_capacity, 2 * buf->m_capacity);
    buf->m_data = realloc(buf->m_data, buf->m_capacity * sizeof(*buf->m_data));
  }
}

static void
particle_buf_push_back(particle_buf_t *buf, particle_t *p)
{
  particle_buf_reserve(buf, buf->m_size + 1);
  buf->m_data[buf->m_size] = *p;
  buf->m_size++;
}

static particle_t *
particle_buf_begin(particle_buf_t *buf)
{
  return buf->m_data;
}

static particle_t *
particle_buf_end(particle_buf_t *buf)
{
  return buf->m_data + buf->m_size;
}

static void
particle_buf_copy(particle_t *from_begin, particle_t *from_end, particle_t *to_begin)
{
  memcpy(to_begin, from_begin, (from_end - from_begin) * sizeof(*to_begin));
}

struct ddcp_info_by_rank {
  struct ddcp_send_entry {
    int patch; // source patch (source rank is this rank)
    int nei_patch; // target patch (target rank is index in send_entry)
    int dir1;  // direction
    int dir1neg;
  } *send_entry;
  int *send_cnts;
  int n_send_entries;
  int n_send;

  struct ddcp_recv_entry { // needs to be same as send_entry with different order!
    int nei_patch;
    int patch;
    int dir1neg;
    int dir1;
  } *recv_entry;
  int *recv_cnts;
  int n_recv_entries;
  int n_recv;

  int rank;
};

struct ddcp_nei {
  particle_buf_t send_buf;
  int n_recv;
  int rank;
  int patch;
};

struct ddcp_patch {
  struct ddcp_nei nei[N_DIR];
  int n_recv;
  ddcp_buf_t buf;
};

struct ddc_particles {
  int nr_patches;
  struct ddcp_patch *patches;
  struct ddcp_info_by_rank *by_rank;
  struct ddcp_info_by_rank *cinfo; // compressed info
  int n_ranks;
  MPI_Request *send_reqs;
  MPI_Request *recv_reqs;

  struct mrc_domain *domain;
};

// ----------------------------------------------------------------------
// ddc_particles_create

static struct ddc_particles *
ddc_particles_create(struct mrc_domain *domain)
{
  struct ddc_particles *ddcp = calloc(1, sizeof(*ddcp));

  ddcp->domain = domain;
  mrc_domain_get_patches(domain, &ddcp->nr_patches);
  ddcp->patches = calloc(ddcp->nr_patches, sizeof(*ddcp->patches));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];

	  particle_buf_ctor(&nei->send_buf);

	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    // use this one as buffer for particles that stay in the same patch
	    nei->rank = -1;
	  } else {
	    mrc_domain_get_neighbor_rank_patch(domain, p, dir, &nei->rank, &nei->patch);
	  }
	}
      }
    }
  }

  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  ddcp->by_rank = calloc(size, sizeof(*ddcp->by_rank));
  struct ddcp_info_by_rank *info = ddcp->by_rank;

  int dir[3];

  // count how many recv_entries per rank
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  info[nei->rank].n_recv_entries++;
	}
      }
    }
  }

  // alloc recv_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      info[r].recv_entry =
	malloc(info[r].n_recv_entries * sizeof(*info[r].recv_entry));
      info[r].recv_cnts =
	malloc(info[r].n_recv_entries * sizeof(*info[r].recv_cnts));
    }
  }

#if 0
  // set up recv_entries
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct recv_entry *re =
	    &info[nei->rank].recv_entry[info[nei->rank].n_recv_entries++];
	  re->patch = p;
	  re->nei_patch = nei->patch;
	  re->dir1 = dir1;
	  re->dir1neg = dir1neg;
	}
      }
    }
  }
#endif

  // count send_entries
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  info[nei->rank].n_send_entries++;
	}
      }
    }
  }

  // alloc send_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      info[r].send_entry =
	malloc(info[r].n_send_entries * sizeof(*info[r].send_entry));
      info[r].send_cnts = 
	malloc(info[r].n_send_entries * sizeof(*info[r].send_cnts));
      info[r].n_send_entries = 0;
    }
  }

  // set up send_entries
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  struct ddcp_send_entry *se =
	    &info[nei->rank].send_entry[info[nei->rank].n_send_entries++];
	  se->patch = p;
	  se->nei_patch = nei->patch;
	  se->dir1 = dir1;
	  se->dir1neg = dir1neg;
	}
      }
    }
  }

  int n_recv_ranks = 0;
  int n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      n_recv_ranks++;
    }
    if (info[r].n_send_entries) {
      n_send_ranks++;
    }
  }

  assert(n_send_ranks == n_recv_ranks);
  ddcp->n_ranks = n_send_ranks;

  ddcp->send_reqs = malloc(ddcp->n_ranks * sizeof(*ddcp->send_reqs));
  ddcp->recv_reqs = malloc(ddcp->n_ranks * sizeof(*ddcp->recv_reqs));

  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      MPI_Irecv(info[r].recv_entry,
		sizeof(struct ddcp_recv_entry) / sizeof(int) * info[r].n_recv_entries,
		MPI_INT, r, 111, comm, &ddcp->recv_reqs[n_recv_ranks++]);
    }
  }  

  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      MPI_Isend(info[r].send_entry,
		sizeof(struct ddcp_send_entry) / sizeof(int) * info[r].n_send_entries,
		MPI_INT, r, 111, comm, &ddcp->send_reqs[n_send_ranks++]);
    }
  }  

  // FIXME / OPT, we're copying alloc'd pointers over,
  // fragile, though correct. info could be free'd here or even
  // earlier
  ddcp->cinfo = malloc(ddcp->n_ranks * sizeof(*ddcp->cinfo));
  struct ddcp_info_by_rank *cinfo = ddcp->cinfo;
  int i = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      assert(info[r].n_send_entries);
      cinfo[i] = info[r];
      cinfo[i].rank = r;
      i++;
    }
  }
  assert(i == ddcp->n_ranks);

  return ddcp;
}

// ----------------------------------------------------------------------
// ddc_particles_destroy

static void
ddc_particles_destroy(struct ddc_particles *ddcp)
{
  if (!ddcp)
    return;

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];

	  particle_buf_dtor(&nei->send_buf);
	}
      }
    }
  }
  free(ddcp->patches);

  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int size;
  MPI_Comm_size(comm, &size);

  struct ddcp_info_by_rank *info = ddcp->by_rank;
  for (int r = 0; r < size; r++) {
    free(info[r].send_entry);
    free(info[r].recv_entry);
    free(info[r].send_cnts);
    free(info[r].recv_cnts);
  }
  free(ddcp->by_rank);
  free(ddcp->cinfo);
  free(ddcp->recv_reqs);
  free(ddcp->send_reqs);

  free(ddcp);
}

// ----------------------------------------------------------------------
// ddc_particles_queue

static void
ddc_particles_queue(struct ddc_particles *ddcp, struct ddcp_patch *patch,
		    int dir[3], particle_t *p)
{
  struct ddcp_nei *nei = &patch->nei[mrc_ddc_dir2idx(dir)];

  particle_buf_push_back(&nei->send_buf, p);
}

// ----------------------------------------------------------------------
// ddc_particles_comm
//
// OPT: could use MPI_Waitany?
// OPT: overall more async
// OPT: 1d instead of 3d loops
// OPT: make the status buffers only as large as needed?

static void
ddc_particles_comm(struct ddc_particles *ddcp, struct psc_mparticles *mprts)
{
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // FIXME, this is assuming our struct is equiv to an array of real_type
  assert(sizeof(particle_t) % sizeof(particle_real_t) == 0);
  int sz = sizeof(particle_t) / sizeof(particle_real_t);
  int dir[3];

  struct ddcp_info_by_rank *cinfo = ddcp->cinfo;

  for (int r = 0; r < ddcp->n_ranks; r++) {
    MPI_Irecv(cinfo[r].recv_cnts, cinfo[r].n_recv_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &ddcp->recv_reqs[r]);
  }  

  for (int r = 0; r < ddcp->n_ranks; r++) {
    cinfo[r].n_send = 0;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      struct ddcp_send_entry *se = &cinfo[r].send_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[se->patch];
      struct ddcp_nei *nei = &patch->nei[se->dir1];
      unsigned int n_send = particle_buf_size(&nei->send_buf);
      cinfo[r].send_cnts[i] = n_send;
      cinfo[r].n_send += n_send;
    }
    MPI_Isend(cinfo[r].send_cnts, cinfo[r].n_send_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &ddcp->send_reqs[r]);
  }  

  // overlap: count local # particles
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->n_recv = 0;
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei_send = &ddcp->patches[nei->patch].nei[dir1neg];
	  patch->n_recv += particle_buf_size(&nei_send->send_buf);
	}
      }
    }
  }

#if 0
  // still want to figure out how to avoid sending info i can calc
  // just need to fix order
  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &info[r].recv_entry[i];
      mprintf("R %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r,
      	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv);
    }
  }

  for (int r = 0; r < size; r++) {
    for (int i = 0; i < info[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &info[r].recv_entry_[i];
      mprintf("r %d %d:%d -> %d:%d dir1 %02d n_recv %d\n", i, r,
      	      re->nei_patch, rank, re->patch, re->dir1, re->n_recv);
    }
  }
#endif

  MPI_Waitall(ddcp->n_ranks, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->n_ranks, ddcp->send_reqs, MPI_STATUSES_IGNORE);

  // add remote # particles
  for (int r = 0; r < ddcp->n_ranks; r++) {
    cinfo[r].n_recv = 0;
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      patch->n_recv += cinfo[r].recv_cnts[i];
      cinfo[r].n_recv += cinfo[r].recv_cnts[i];
    }
  }

  // realloc
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    ddcp_buf_reserve(&patch->buf, ddcp_buf_size(&patch->buf) + patch->n_recv);
  }

  // leave room for receives (FIXME? just change order)
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      ddcp_buf_resize(&patch->buf, ddcp_buf_size(&patch->buf) + cinfo[r].recv_cnts[i]);
    }
  }

  int n_send = 0, n_recv = 0;
  for (int r = 0; r < ddcp->n_ranks; r++) {
    n_send += cinfo[r].n_send;
    n_recv += cinfo[r].n_recv;
  }
  
  // post sends
  particle_buf_t send_buf;
  particle_buf_ctor(&send_buf);
  particle_buf_reserve(&send_buf, n_send);
  particle_t *it = particle_buf_begin(&send_buf);
  for (int r = 0; r < ddcp->n_ranks; r++) {
    if (cinfo[r].n_send == 0)
      continue;

    particle_t *it0 = it;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      struct ddcp_send_entry *se = &cinfo[r].send_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[se->patch];
      particle_buf_t *send_buf_nei = &patch->nei[se->dir1].send_buf;
      particle_buf_copy(particle_buf_begin(send_buf_nei), particle_buf_end(send_buf_nei), it);
      it += particle_buf_size(send_buf_nei);
    }
    MPI_Isend(it0, sz * cinfo[r].n_send, MPI_PARTICLES_REAL,
	      cinfo[r].rank, 1, comm, &ddcp->send_reqs[r]);
  }
  assert(it == particle_buf_begin(&send_buf) + n_send);

  // post receives
  particle_buf_t recv_buf;
  particle_buf_ctor(&recv_buf);
  particle_buf_reserve(&recv_buf, n_recv);
  it = particle_buf_begin(&recv_buf);
  for (int r = 0; r < ddcp->n_ranks; r++) {
    if (cinfo[r].n_recv == 0)
      continue;

    MPI_Irecv(it, sz * cinfo[r].n_recv, MPI_PARTICLES_REAL,
	      cinfo[r].rank, 1, comm, &ddcp->recv_reqs[r]);
    it += cinfo[r].n_recv;
  }
  assert(it == particle_buf_begin(&recv_buf) + n_recv);

  // overlap: copy particles from local proc
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dir1neg = mrc_ddc_dir2idx((int[3]) { -dir[0], -dir[1], -dir[2] });
	  struct ddcp_nei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  particle_buf_t *nei_send_buf = &ddcp->patches[nei->patch].nei[dir1neg].send_buf;
	  unsigned int end = ddcp_buf_size(&patch->buf);
	  ddcp_buf_resize(&patch->buf, end + particle_buf_size(nei_send_buf));
	  particle_buf_copy(particle_buf_begin(nei_send_buf), particle_buf_end(nei_send_buf),
			    ddcp_buf_at(&patch->buf, end));
	}
      }
    }
  }

  MPI_Waitall(ddcp->n_ranks, ddcp->recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ddcp->n_ranks, ddcp->send_reqs, MPI_STATUSES_IGNORE);

  // copy received particles into right place

  int *old_head = malloc(ddcp->nr_patches * sizeof(*old_head));
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    old_head[p] = ddcp_buf_size(&patch->buf);
    ddcp_buf_resize(&patch->buf, old_head[p] - patch->n_recv);
  }

  it = particle_buf_begin(&recv_buf);
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      int end = ddcp_buf_size(&patch->buf);
      ddcp_buf_resize(&patch->buf, end + cinfo[r].recv_cnts[i]);
      particle_buf_copy(it, it + cinfo[r].recv_cnts[i], ddcp_buf_at(&patch->buf, end));
      it += cinfo[r].recv_cnts[i];
    }
  }
  assert(it == particle_buf_begin(&recv_buf) + n_recv);

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    ddcp_buf_resize(&patch->buf, old_head[p]);
  }
  free(old_head);

  particle_buf_dtor(&send_buf);
  particle_buf_dtor(&recv_buf);

  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp_buf_dtor(&ddcp->patches[p].buf);
  }
}


// ----------------------------------------------------------------------
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  bnd->ddcp = ddc_particles_create(bnd->psc->mrc_domain);

#if DDCP_TYPE == DDCP_TYPE_COMMON
  psc_bnd_particles_open_setup(bnd);
#endif
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);

#if DDCP_TYPE == DDCP_TYPE_COMMON
  psc_bnd_particles_open_unsetup(bnd);
#endif
}

#if DDCP_TYPE == DDCP_TYPE_COMMON2

// ----------------------------------------------------------------------
// find_block_indices_count

static inline void
find_block_indices_count(unsigned int *b_idx, unsigned int *b_cnts,
			 struct psc_mparticles *mprts, int p, int off)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];

  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  int *b_mx = patch->b_mx;
  for (int i = off; i < n_prts; i++) {
    particle_t *part = particle_iter_at(prts.begin, i);
    int b_pos[3];
    particle_xi_get_block_pos(&part->xi, patch->b_dxi, b_pos);
    assert(b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	   b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	   b_pos[2] >= 0 && b_pos[2] < b_mx[2]);
    b_idx[i] = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    b_cnts[b_idx[i]]++;
  }
}

// ----------------------------------------------------------------------
// find_block_indices_count_reorder

static void _mrc_unused
find_block_indices_count_reorder(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  unsigned int cnt = n_prts;
  int *b_mx = patch->b_mx;
  memset(patch->b_cnt, 0, (patch->nr_blocks + 1) * sizeof(*patch->b_cnt));

  for (int i = 0; i < n_prts; i++) {
    particle_t *part = particle_iter_at(prts.begin, i);
    int b_pos[3];
    particle_xi_get_block_pos(&part->xi, patch->b_dxi, b_pos);
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      patch->b_idx[i] = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    } else { // out of bounds
      patch->b_idx[i] = patch->nr_blocks;
      *particle_iter_at(prts.begin, cnt) = *part;
      cnt++;
    }
    patch->b_cnt[patch->b_idx[i]]++;
  }
}

static void _mrc_unused
count_and_reorder_to_back(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  particle_range_t prts = particle_range_mprts(mprts, p);

  memset(patch->b_cnt, 0, (patch->nr_blocks + 1) * sizeof(*patch->b_cnt));
  unsigned int n_prts = particle_range_size(prts);
  unsigned int cnt = n_prts;
  for (int i = 0; i < n_prts; i++) {
    if (patch->b_idx[i] == patch->nr_blocks) {
      *particle_iter_at(prts.begin, cnt) = *particle_iter_at(prts.begin, i);
      cnt++;
    }
    patch->b_cnt[patch->b_idx[i]]++;
  }
}

static void _mrc_unused
reorder_to_back(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  particle_range_t prts = particle_range_mprts(mprts, p);

  unsigned int n_prts = particle_range_size(prts);
  unsigned int cnt = n_prts;
  for (int i = 0; i < n_prts; i++) {
    if (patch->b_idx[i] == patch->nr_blocks) {
      *particle_iter_at(prts.begin, cnt) = *particle_iter_at(prts.begin, i);
      cnt++;
    }
  }
}

// ----------------------------------------------------------------------
// count_block_indices

static inline void
count_block_indices(unsigned int *b_cnts, unsigned int *b_idx, int n_prts, int off)
{
  for (int i = off; i < n_prts; i++) {
    b_cnts[b_idx[i]]++;
  }
}

// ----------------------------------------------------------------------
// exclusive_scan

static inline void
exclusive_scan(unsigned int *b_cnts, int n)
{
  unsigned int sum = 0;
  for (int i = 0; i < n; i++) {
    unsigned int cnt = b_cnts[i];
    b_cnts[i] = sum;
    sum += cnt;
  }
}

// ----------------------------------------------------------------------
// sort_indices

static void
sort_indices(unsigned int *b_idx, unsigned int *b_sum, unsigned int *b_ids, int n_prts)
{
  for (int n = 0; n < n_prts; n++) {
    unsigned int n_new = b_sum[b_idx[n]]++;
    assert(n_new < n_prts);
    b_ids[n_new] = n;
  }
}

#endif

// ======================================================================

#if DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_CUDA

static inline int
get_n_send(struct psc_mparticles *mprts, int p)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON2
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  return sub->patch[p].n_send;
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  return cmprts->bnd.bpatch[p].n_send;
#endif
}

static inline int
get_head(struct psc_mparticles *mprts, int p)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON2
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  particle_range_t prts = particle_range_mprts(mprts, p);
  return particle_range_size(prts) - sub->patch[p].n_send;
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  return 0;
#endif
}

static void
exchange_particles_pre(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts, int p)
{
  struct psc *psc = bnd->psc;
  struct ddc_particles *ddcp = bnd->ddcp;

  int n_send = get_n_send(mprts, p);
  struct psc_patch *patch = &psc->patch[p];
  particle_real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = patch->ldims[d] * patch->dx[d];
  }
  const particle_real_t *b_dxi = mparticles_patch_get_b_dxi(mprts, p);
  const int *b_mx = mparticles_patch_get_b_mx(mprts, p);
  
  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  struct ddcp_patch *ddcp_patch = &ddcp->patches[p];
  ddcp_buf_ctor(&ddcp_patch->buf, mprts, p);
  ddcp_buf_resize(&ddcp_patch->buf, get_head(mprts, p));
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    particle_buf_resize(&ddcp_patch->nei[dir1].send_buf, 0);
  }

  int n_end = ddcp_buf_size(&ddcp_patch->buf) + n_send;
  for (int n = ddcp_buf_size(&ddcp_patch->buf); n < n_end; n++) {
    particle_t *prt = ddcp_buf_at(&ddcp_patch->buf, n);
    particle_real_t *xi = &prt->xi;
    particle_real_t *pxi = &prt->pxi;
    
    bool drop = false;
    int dir[3];
    for (int d = 0; d < 3; d++) {
      int bi = particle_real_fint(xi[d] * b_dxi[d]);
      if (bi < 0) {
	// FIXME, assumes every patch has same dimensions
	if (patch->off[d] != 0 || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	  xi[d] += xm[d];
	  dir[d] = -1;
	  bi = particle_real_fint(xi[d] * b_dxi[d]);
	  if (bi >= b_mx[d]) {
	    xi[d] = 0.;
	    dir[d] = 0;
	  }
	} else {
	  switch (psc->domain.bnd_part_lo[d]) {
	  case BND_PART_REFLECTING:
	    xi[d]  = -xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    break;
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else if (bi >= b_mx[d]) {
	if (patch->off[d] + patch->ldims[d] != psc->domain.gdims[d] ||
	    psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	  xi[d] -= xm[d];
	  dir[d] = +1;
	  bi = particle_real_fint(xi[d] * b_dxi[d]);
	  if (bi < 0) {
	    xi[d] = 0.;
	  }
	} else {
	  switch (psc->domain.bnd_part_hi[d]) {
	  case BND_PART_REFLECTING:
	    xi[d] = 2.f * xm[d] - xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] *= (1. - 1e-6);
	    }
	    break;
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else {
	dir[d] = 0;
      }
      if (!drop) {
	if (xi[d] < 0.f && xi[d] > -1e-6f) {
	  mprintf("d %d xi %g\n", d, xi[d]);
	  xi[d] = 0.f;
	}
	assert(xi[d] >= 0.f);
	assert(xi[d] <= xm[d]);
      }
    }
    if (!drop) {
      if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	ddcp_buf_push_back(&ddcp_patch->buf, prt);
      } else {
	ddc_particles_queue(ddcp, ddcp_patch, dir, prt);
      }
    }
  }
}

#endif

#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP || DDCP_TYPE == DDCP_TYPE_COMMON2

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_prep

static void
psc_bnd_particles_sub_exchange_particles_prep(struct psc_bnd_particles *bnd,
					      struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = particle_range_mprts(mprts, p);

#if DDCP_TYPE == DDCP_TYPE_COMMON2
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  if (1) {
    //      find_block_indices_count_reorderx(prts);
    count_and_reorder_to_back(mprts, p);
  }
  patch->n_part_save = particle_range_size(prts);
  patch->n_send = patch->b_cnt[patch->nr_blocks];
  mparticles_patch_resize(mprts, p, mparticles_get_n_prts(mprts, p) + patch->n_send);

  exchange_particles_pre(bnd, mprts, p);

#elif DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP

  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc *psc = bnd->psc;

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  struct psc_patch *ppatch = &psc->patch[p];
  particle_real_t b_dxi[3] = { 1.f / ppatch->dx[0], 1.f / ppatch->dx[1], 1.f / ppatch->dx[2] };
  particle_real_t xm[3];
  int b_mx[3];
  for (int d = 0; d < 3; d++ ) {
    b_mx[d] = ppatch->ldims[d];
    xm[d] = b_mx[d] * ppatch->dx[d];
  }
  
  struct ddcp_patch *patch = &ddcp->patches[p];
  ddcp_buf_ctor(&patch->buf, mprts, p);
  ddcp_buf_resize(&patch->buf, 0);
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    particle_buf_resize(&patch->nei[dir1].send_buf, 0);
  }
  unsigned int n_prts = particle_range_size(prts);
  for (int i = 0; i < n_prts; i++) {
    particle_t *part = particle_iter_at(prts.begin, i);
    particle_real_t *xi = &part->xi; // slightly hacky relies on xi, yi, zi to be contiguous in the struct. FIXME
    particle_real_t *pxi = &part->pxi;
    
    int b_pos[3];
    particle_xi_get_block_pos(xi, b_dxi, b_pos);
    if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] && // OPT, could be optimized with casts to unsigned
	b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      // fast path
      // inside domain: move into right position
      ddcp_buf_push_back(&patch->buf, part);
    } else {
      // slow path
      bool drop = false;
      int dir[3];
      for (int d = 0; d < 3; d++) {
	if (b_pos[d] < 0) {
	  if (!at_lo_boundary(p, d) || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	    xi[d] += xm[d];
	    dir[d] = -1;
	    int bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] = 0.;
	      dir[d] = 0;
	    }
	  } else {
	    switch (psc->domain.bnd_part_lo[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] =  -xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else if (b_pos[d] >= b_mx[d]) {
	  if (!at_hi_boundary(p, d) ||
	      psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	    xi[d] -= xm[d];
	    dir[d] = +1;
	    int bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi < 0) {
	      xi[d] = 0.;
	    }
	  } else {
	    switch (psc->domain.bnd_part_hi[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xm[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      int bi = particle_real_fint(xi[d] * b_dxi[d]);
	      if (bi >= b_mx[d]) {
		xi[d] *= (1. - 1e-6);
	      }
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else {
	  // computational bnd
	  dir[d] = 0;
	}
	if (!drop) {
	  if (xi[d] < 0.f && xi[d] > -1e-6f) {
	    //	    mprintf("d %d xi %g\n", d, xi[d]);
	    xi[d] = 0.f;
	  }
	  assert(xi[d] >= 0.f);
	  assert(xi[d] <= xm[d]);
	}
      }
      if (!drop) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  ddcp_buf_push_back(&patch->buf, part);
	} else {
	  ddc_particles_queue(ddcp, patch, dir, part);
	}
      }

    }
  }
#endif
}

#endif

#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP || DDCP_TYPE == DDCP_TYPE_COMMON2

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_post

static void
psc_bnd_particles_sub_exchange_particles_post(struct psc_bnd_particles *bnd,
					      struct psc_mparticles *mprts, int p)
{
  struct ddc_particles *ddcp = bnd->ddcp;
  struct ddcp_patch *dpatch = &ddcp->patches[p];

  mparticles_patch_resize(mprts, p, ddcp_buf_size(&dpatch->buf));
  
#if DDCP_TYPE == DDCP_TYPE_COMMON2
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  particle_range_t prts = particle_range_mprts(mprts, p);

  find_block_indices_count(patch->b_idx, patch->b_cnt, mprts, p, patch->n_part_save);
  exclusive_scan(patch->b_cnt, patch->nr_blocks + 1);
  sort_indices(patch->b_idx, patch->b_cnt, patch->b_ids, particle_range_size(prts));
  
  mparticles_patch_resize(mprts, p, patch->b_cnt[patch->nr_blocks - 1]);
  patch->need_reorder = true; // FIXME, need to honor before get()/put()
#endif
}

#endif

#if DDCP_TYPE == DDCP_TYPE_CUDA

// ----------------------------------------------------------------------
// mprts_convert_to_cuda

static void
mprts_convert_to_cuda(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int nr_recv = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct ddc_particles *ddcp = bnd->ddcp;
    struct ddcp_patch *patch = &ddcp->patches[p];
    nr_recv += ddcp_buf_size(&patch->buf);
  }

  cmprts->bnd.h_bnd_xi4  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_xi4));
  cmprts->bnd.h_bnd_pxi4 = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_pxi4));
  cmprts->bnd.h_bnd_idx  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_idx));
  cmprts->bnd.h_bnd_off  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_off));

  memset(cmprts->bnd.h_bnd_cnt, 0,
	 cmprts->n_blocks * sizeof(*cmprts->bnd.h_bnd_cnt));

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct ddc_particles *ddcp = bnd->ddcp;
    struct ddcp_patch *patch = &ddcp->patches[p];
    int n_recv = ddcp_buf_size(&patch->buf);
    cmprts->bnd.bpatch[p].n_recv = n_recv;
    
    float4 *h_bnd_xi4 = cmprts->bnd.h_bnd_xi4 + off;
    float4 *h_bnd_pxi4 = cmprts->bnd.h_bnd_pxi4 + off;
    unsigned int *h_bnd_idx = cmprts->bnd.h_bnd_idx + off;
    unsigned int *h_bnd_off = cmprts->bnd.h_bnd_off + off;
    for (int n = 0; n < n_recv; n++) {
      particle_t *prt = &cmprts->bnd.bpatch[p].prts[n];
      h_bnd_xi4[n].x  = prt->xi;
      h_bnd_xi4[n].y  = prt->yi;
      h_bnd_xi4[n].z  = prt->zi;
      h_bnd_xi4[n].w  = cuda_int_as_float(prt->kind);
      h_bnd_pxi4[n].x = prt->pxi;
      h_bnd_pxi4[n].y = prt->pyi;
      h_bnd_pxi4[n].z = prt->pzi;
      h_bnd_pxi4[n].w = prt->qni_wni;

      int b_pos[3];
      for (int d = 0; d < 3; d++) {
	float *xi = &h_bnd_xi4[n].x;
	b_pos[d] = particle_real_fint(xi[d] * cmprts->b_dxi[d]);
	if (b_pos[d] < 0 || b_pos[d] >= cmprts->b_mx[d]) {
	  printf("!!! xi %g %g %g\n", xi[0], xi[1], xi[2]);
	  printf("!!! d %d xi4[n] %g biy %d // %d\n",
		 d, xi[d], b_pos[d], cmprts->b_mx[d]);
	  if (b_pos[d] < 0) {
	    xi[d] = 0.f;
	  } else {
	    xi[d] *= (1. - 1e-6);
	  }
	}
	b_pos[d] = particle_real_fint(xi[d] * cmprts->b_dxi[d]);
	assert(b_pos[d] >= 0 && b_pos[d] < cmprts->b_mx[d]);
      }
      unsigned int b = (b_pos[2] * cmprts->b_mx[1] + b_pos[1]) * cmprts->b_mx[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n] = b;
      h_bnd_off[n] = cmprts->bnd.h_bnd_cnt[b]++;
    }
    free(cmprts->bnd.bpatch[p].prts);
    cmprts->bnd.bpatch[p].capacity = 0;
    off += n_recv;
  }
}
#endif

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_prep

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

static void
psc_bnd_particles_sub_exchange_mprts_prep(struct psc_bnd_particles *bnd,
					  struct psc_mparticles *mprts)
{
  prof_restart(pr_time_step_no_comm);

  static int pr_A, pr_B, pr_D, pr_E, pr_F, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0= prof_register("xchg_reduce", 1., 0, 0);
    pr_B1= prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
    pr_E = prof_register("xchg_cvt_from", 1., 0, 0);
    pr_F = prof_register("xchg_pre", 1., 0, 0);
  }

#if DDCP_TYPE == DDCP_TYPE_CUDA

  //prof_start(pr_A);
  //cuda_mprts_find_block_keys(mprts);
  //prof_stop(pr_A);
  
  prof_start(pr_B0);
  cuda_mprts_spine_reduce(mprts);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  cuda_mprts_find_n_send(mprts);
  prof_stop(pr_B1);

  prof_start(pr_B);
  cuda_mprts_scan_send_buf_total(mprts);
  prof_stop(pr_B);

  prof_start(pr_D);
  cuda_mprts_copy_from_dev(mprts);
  prof_stop(pr_D);
  
  prof_start(pr_E);
  cuda_mprts_convert_from_cuda(mprts);
  prof_stop(pr_E);
#endif
  
  prof_start(pr_F);
#if DDCP_TYPE == DDCP_TYPE_COMMON
#pragma omp parallel for
#endif
  for (int p = 0; p < mprts->nr_patches; p++) {
    psc_balance_comp_time_by_patch[p] -= MPI_Wtime();
#if DDCP_TYPE == DDCP_TYPE_CUDA
    exchange_particles_pre(bnd, mprts, p);
#elif DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP || DDCP_TYPE == DDCP_TYPE_COMMON2
    psc_bnd_particles_sub_exchange_particles_prep(bnd, mprts, p);
#endif
    psc_balance_comp_time_by_patch[p] += MPI_Wtime();
  }
  prof_stop(pr_F);

  prof_stop(pr_time_step_no_comm);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_post

static void
psc_bnd_particles_sub_exchange_mprts_post(struct psc_bnd_particles *bnd,
					  struct psc_mparticles *mprts)
{
#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP || DDCP_TYPE == DDCP_TYPE_COMMON2
  for (int p = 0; p < mprts->nr_patches; p++) {
    psc_bnd_particles_sub_exchange_particles_post(bnd, mprts, p);
  }
#elif DDCP_TYPE == DDCP_TYPE_CUDA
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  static int pr_A, pr_B, pr_C, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_cvt_to", 1., 0, 0);
    pr_B = prof_register("xchg_to_dev", 1., 0, 0);
    pr_C = prof_register("xchg_bidx", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1= prof_register("xchg_upd_off", 1., 0, 0);
    pr_E = prof_register("xchg_reorder", 1., 0, 0);
  }

  prof_start(pr_A);
  mprts_convert_to_cuda(bnd, mprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  cuda_mprts_copy_to_dev(mprts);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mprts_find_block_indices_3(mprts);
  prof_stop(pr_C);
  
  prof_start(pr_D);
  unsigned int n_prts_by_patch[mprts->nr_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  cuda_mprts_sort(mprts, (int *) n_prts_by_patch); // FIXME cast
  // FIXME, is this necessary, or doesn't update_offsets() do this, too?
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  prof_stop(pr_D);

  prof_start(pr_D1);
  cuda_mprts_update_offsets(mprts);
  prof_stop(pr_D1);
  
#if 0
  prof_start(pr_E);
  cuda_mparticles_reorder(cmprts);
  prof_stop(pr_E);
  //  cuda_mprts_check_ordered_total(mprts);
#else
  cmprts->need_reorder = true;
#endif
#endif
}

#if DDCP_TYPE == DDCP_TYPE_CUDA
// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_serial_periodic
//
// specialized version if there's only one single patch,
// with all periodic b.c.
//
// TODO: a lot of consolidation with the generic case should
// be possible. In particular, optimizations carry over, like
// calculating new offsets as part of the spine, not calculating
// new block indices, not calculates old ids.
// The most significant stumbling block to make them pretty much the same
// is that this case here needs to handle the odd periodic spine.

static void
psc_bnd_particles_sub_exchange_particles_serial_periodic(struct psc_bnd_particles *psc_bnd_particles,
						struct psc_mparticles *mprts)
{
  assert(0);
#if 0
  static int pr_F, pr_G, pr_H;
  if (!pr_F) {
    pr_F = prof_register("xchg_bidx_ids", 1., 0, 0);
    pr_G = prof_register("xchg_sort_pairs", 1., 0, 0);
    pr_H = prof_register("xchg_reorder_off", 1., 0, 0);
  }

  cuda_exchange_particles(0, psc_mparticles_get_patch(particles, 0));

  // sort
  for (int p = 0; p < particles->nr_patches; p++) {
    prof_start(pr_F);
    cuda_find_block_indices(prts, cuda->h_dev->bidx);
    prof_stop(pr_F);

    prof_start(pr_G);
    sort_pairs_device_2(cuda->sort_ctx, cuda->h_dev->bidx,
			cuda->h_dev->alt_ids,
			mprts_cuda->n_prts_by_patch[p],
			cuda->h_dev->offsets);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->h_dev->alt_ids);
    prof_stop(pr_H);
  }
#endif
}

#endif

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_general

static void
psc_bnd_particles_sub_exchange_particles_general(struct psc_bnd_particles *bnd,
						 struct psc_mparticles *mprts)
{
  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep", 1., 0, 0);
    pr_B = prof_register("xchg_comm", 1., 0, 0);
    pr_C = prof_register("xchg_post", 1., 0, 0);
  }
  
  struct ddc_particles *ddcp = bnd->ddcp;

  prof_start(pr_A);
  psc_bnd_particles_sub_exchange_mprts_prep(bnd, mprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, mprts);
  prof_stop(pr_B);
  
  prof_start(pr_C);
  psc_bnd_particles_sub_exchange_mprts_post(bnd, mprts);
  prof_stop(pr_C);

  //struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", JXI, JXI + 3);
  //psc_bnd_particles_open_boundary(bnd, particles, mflds);
  //psc_mfields_put_as(mflds, psc->flds, JXI, JXI + 3);
}


// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd,
			       struct psc_mparticles *mprts_base)
{
#if DDCP_TYPE == DDCP_TYPE_CUDA
  // This function only makes sense if it's called for particles already being of cuda
  // type. If particles aren't in the right patches, the conversion in get_as would fail...

  assert(strcmp(psc_mparticles_type(mprts_base), "cuda") == 0);
#endif
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  
#if DDCP_TYPE == DDCP_TYPE_CUDA
  int size;
  MPI_Comm_size(psc_bnd_particles_comm(bnd), &size);

  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    psc_bnd_particles_sub_exchange_particles_serial_periodic(bnd, mprts);
  } else {
    psc_bnd_particles_sub_exchange_particles_general(bnd, mprts);
  }
#elif DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON_OMP || DDCP_TYPE == DDCP_TYPE_COMMON2
  psc_bnd_particles_sub_exchange_particles_general(bnd, mprts);
#endif

  psc_mparticles_put_as(mprts, mprts_base, 0);
}
