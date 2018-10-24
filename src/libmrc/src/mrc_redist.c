
#include <mrc_redist.h>
#include <mrc_domain.h>

#include <stdlib.h>

void
mrc_redist_init(struct mrc_redist *redist, struct mrc_domain *domain,
		int slab_offs[3], int slab_dims[3], int nr_writers)
{
  redist->domain = domain;
  redist->comm = mrc_domain_comm(domain);
  MPI_Comm_rank(redist->comm, &redist->rank);
  MPI_Comm_size(redist->comm, &redist->size);

  redist->nr_writers = nr_writers;
  redist->writer_ranks = calloc(nr_writers, sizeof(*redist->writer_ranks));
  // setup writers, just use first nr_writers ranks,
  // could do something fancier in the future
  redist->is_writer = 0;
  for (int i = 0; i < nr_writers; i++) {
    redist->writer_ranks[i] = i;
    if (i == redist->rank) {
      redist->is_writer = 1;
    }
  }
  
  MPI_Comm_split(redist->comm, redist->is_writer, redist->rank, &redist->comm_writers);

  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);
  int nr_global_patches;
  mrc_domain_get_nr_global_patches(domain, &nr_global_patches);
  for (int d = 0; d < 3; d++) {
    if (slab_dims[d]) {
      redist->slab_dims[d] = slab_dims[d];
    } else {
      redist->slab_dims[d] = gdims[d];
    }
    redist->slab_offs[d] = slab_offs[d];
  }
  redist->slow_dim = 2;
  while (gdims[redist->slow_dim] == 1) {
    redist->slow_dim--;
  }
  assert(redist->slow_dim >= 0);
  int total_slow_indices = redist->slab_dims[redist->slow_dim];
  redist->slow_indices_per_writer = total_slow_indices / nr_writers;
  redist->slow_indices_rmndr = total_slow_indices % nr_writers;
}

void
mrc_redist_destroy(struct mrc_redist *redist)
{
  free(redist->writer_ranks);
  MPI_Comm_free(&redist->comm_writers);
}

static void
mrc_redist_writer_offs_dims(struct mrc_redist *redist, int writer,
			    int *writer_offs, int *writer_dims)
{
  for (int d = 0; d < 3; d++) {
    writer_dims[d] = redist->slab_dims[d];
    writer_offs[d] = redist->slab_offs[d];
  }
  writer_dims[redist->slow_dim] = redist->slow_indices_per_writer + (writer < redist->slow_indices_rmndr);
  if (writer < redist->slow_indices_rmndr) {
    writer_offs[redist->slow_dim] += (redist->slow_indices_per_writer + 1) * writer;
  } else {
    writer_offs[redist->slow_dim] += redist->slow_indices_rmndr +
      redist->slow_indices_per_writer * writer;
  }
}

#define BUFLOOP(ix, iy, iz, ilo, hi) \
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {\
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {\
	  for (int ix = ilo[0]; ix < ihi[0]; ix++)

#define BUFLOOP_END }}

MPI_Datatype to_mpi_datatype(int data_type)
{
  switch (data_type) {
  case MRC_NT_FLOAT:  return MPI_FLOAT; break;
  case MRC_NT_DOUBLE: return MPI_DOUBLE; break;
  case MRC_NT_INT:    return MPI_INT; break;
  default: assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_redist_write_send_init

static void
mrc_redist_write_send_init(struct mrc_redist *redist, int size_of_type)
{
  struct mrc_redist_write_send *send = &redist->write_send;

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(redist->domain, &nr_patches);

  // count number of writers we actually need to communicate with
  int n_peers = 0;
  for (int writer = 0; writer < redist->nr_writers; writer++) {
    // don't send to self
    int writer_rank = redist->writer_ranks[writer];
    if (writer_rank == redist->rank) {
      continue;
    }
    
    int writer_offs[3], writer_dims[3];
    mrc_redist_writer_offs_dims(redist, writer, writer_offs, writer_dims);

    // count blocks that we need to send to writer
    int n_blocks = 0;
    for (int p = 0; p < nr_patches; p++) {
      int ilo[3], ihi[3];
      bool has_intersection = find_intersection(ilo, ihi, patches[p].off, patches[p].ldims,
						writer_offs, writer_dims);
      if (!has_intersection) {
	continue;
      }
      n_blocks++;
    }

    if (!n_blocks) {
      continue;
    }

    n_peers++;
  }

  send->buf_size = 0;
  send->peers_begin = calloc(n_peers, sizeof(*send->peers_begin));
  send->peers_end = send->peers_begin + n_peers;
  send->reqs = calloc(n_peers, sizeof(*send->reqs));
  send->disps = calloc(redist->size, sizeof(*send->disps));
  send->cnts = calloc(redist->size, sizeof(*send->cnts));

  struct mrc_redist_peer *w = send->peers_begin;
  for (int writer = 0; writer < redist->nr_writers; writer++) {
    // don't send to self
    int writer_rank = redist->writer_ranks[writer];
    if (writer_rank == redist->rank) {
      continue;
    }
    
    int writer_offs[3], writer_dims[3];
    mrc_redist_writer_offs_dims(redist, writer, writer_offs, writer_dims);

    // count blocks that we need to send to writer
    int n_blocks = 0;
    for (int p = 0; p < nr_patches; p++) {
      int ilo[3], ihi[3];
      bool has_intersection = find_intersection(ilo, ihi, patches[p].off, patches[p].ldims,
						writer_offs, writer_dims);
      if (!has_intersection) {
	continue;
      }
      n_blocks++;
    }

    if (!n_blocks) {
      continue;
    }

    w->rank = writer_rank;
    w->blocks_begin = calloc(n_blocks, sizeof(*w->blocks_begin));
    w->blocks_end = w->blocks_begin + n_blocks;

    w->buf_size = 0;
    struct mrc_redist_block *b = w->blocks_begin;
    for (int p = 0; p < nr_patches; p++) {
      int ilo[3], ihi[3];
      bool has_intersection = find_intersection(ilo, ihi, patches[p].off, patches[p].ldims,
						writer_offs, writer_dims);
      if (!has_intersection) {
	continue;
      }

      b->p = p;
      for (int d = 0; d < 3; d++) {
	b->ilo[d] = ilo[d];
	b->ihi[d] = ihi[d];
      }

      b++;
      w->buf_size += (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);
    }

    // allocate buf per writer
    //mprintf("to writer %d buf_size %d n_blocks %d\n", writer, w->buf_size, w->n_blocks);
    send->cnts[w->rank] = w->buf_size;
    send->buf_size += w->buf_size;
    w++;
  }

  int last = 0;
  for (int r = 0; r < redist->size; r++) {
    send->disps[r] = last;
    last += send->cnts[r];
  }
  assert(last == send->buf_size);
  
  send->buf = malloc(send->buf_size * size_of_type);
  assert(send->buf);
  int off = 0;
  for (struct mrc_redist_peer* w = send->peers_begin; w != send->peers_end; w++) {
    assert(off == send->disps[w->rank]);
    w->off = off;
    off += w->buf_size;
  }

#if 0
  size_t g_data[2], data[2] = { send->buf_size, send->peers_end - send->peers_begin };
  MPI_Reduce(data, g_data, 2, MPI_LONG, MPI_SUM, 0, redist->comm);
  mpi_printf(redist->comm, "avg send buf size = %ld (n peers %ld)\n", g_data[0] / redist->size,
	     g_data[1] / redist->size);
#endif
}

// ----------------------------------------------------------------------
// mrc_redist_write_send_destroy

static void
mrc_redist_write_send_destroy(struct mrc_redist *redist)
{
  struct mrc_redist_write_send *send = &redist->write_send;

  for (struct mrc_redist_peer *w = send->peers_begin; w != send->peers_end; w++) {
    free(w->blocks_begin);
  }
  free(send->peers_begin);
  free(send->reqs);
  free(send->disps);
  free(send->cnts);
  free(send->buf);
}

// ----------------------------------------------------------------------
// mrc_redist_write_send_prep

static void
mrc_redist_write_send_prep(struct mrc_redist *redist, struct mrc_fld *m3, int m)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, &nr_patches);

  struct mrc_redist_write_send *send = &redist->write_send;

  for (struct mrc_redist_peer *w = send->peers_begin; w != send->peers_end; w++) {
    // fill buf per writer
    void *buf = send->buf + w->off * m3->_nd->size_of_type;
    for (struct mrc_redist_block *b = w->blocks_begin; b != w->blocks_end; b++) {
      int p = b->p;
      int *off = patches[p].off;
      int *ilo = b->ilo, *ihi = b->ihi;
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->_domain, p, &info);
      assert(!m3->_aos);
      switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT:
      {
      	float *buf_ptr = (float *) buf;
      	BUFLOOP(ix, iy, iz, ilo, ihi) {
    	    *buf_ptr++ = MRC_S5(m3, ix-off[0],iy-off[1],iz-off[2], m, p);
      	} BUFLOOP_END;
	buf = buf_ptr;
      	break;
      }
      case MRC_NT_DOUBLE:
      {
      	double *buf_ptr = (double *) buf;
      	BUFLOOP(ix, iy, iz, ilo, ihi) {
          *buf_ptr++ = MRC_D5(m3, ix-off[0],iy-off[1],iz-off[2], m, p);
      	} BUFLOOP_END;
	buf = buf_ptr;
      	break;
      }
      case MRC_NT_INT:
      {
      	int *buf_ptr = (int *) buf;
      	BUFLOOP(ix, iy, iz, ilo, ihi) {
    	    *buf_ptr++ = MRC_I5(m3, ix-off[0],iy-off[1],iz-off[2], m, p);
      	} BUFLOOP_END;
	buf = buf_ptr;
      	break;
      }
      default:
      {
      	assert(0);
      }
      }
    }
  }
}

// ----------------------------------------------------------------------
// mrc_redist_write_recv_init

static void
mrc_redist_write_recv_init(struct mrc_redist *redist, int size_of_type)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;

  recv->disps = calloc(redist->size, sizeof(*recv->disps));
  recv->cnts = calloc(redist->size, sizeof(*recv->cnts));
}

// ----------------------------------------------------------------------
// mrc_redist_write_recv_init2

static void
mrc_redist_write_recv_init2(struct mrc_redist *redist, struct mrc_ndarray *nd,
			    int size_of_type)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;

  // find out who's sending, OPT: this way is not very scalable
  // could also be optimized by just looking at slow_dim
  // FIXME, figure out pattern and cache, at least across components

  int nr_global_patches;
  mrc_domain_get_nr_global_patches(redist->domain, &nr_global_patches);

  recv->n_recv_patches = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(redist->domain, gp, &info);

    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     mrc_ndarray_offs(nd), mrc_ndarray_dims(nd));
    if (!has_intersection) {
      continue;
    }

    recv->n_recv_patches++;
  }
  //mprintf("n_recv_patches %d\n", recv->n_recv_patches);

  recv->recv_patches = calloc(recv->n_recv_patches, sizeof(*recv->recv_patches));

  struct mrc_redist_block **recv_patches_by_rank = calloc(redist->size + 1, sizeof(*recv_patches_by_rank));

  int cur_rank = -1;
  recv->n_recv_patches = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(redist->domain, gp, &info);

    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     mrc_ndarray_offs(nd), mrc_ndarray_dims(nd));
    if (!has_intersection) {
      continue;
    }

    assert(info.rank >= cur_rank);
    while (cur_rank < info.rank) {
      cur_rank++;
      recv_patches_by_rank[cur_rank] = &recv->recv_patches[recv->n_recv_patches];
      //mprintf("rank %d patches start at %d\n", cur_rank, recv->n_recv_patches);
    }

    struct mrc_redist_block *recv_patch = &recv->recv_patches[recv->n_recv_patches];
    for (int d = 0; d < 3; d++) {
      recv_patch->ilo[d] = ilo[d];
      recv_patch->ihi[d] = ihi[d];
    }
    recv->n_recv_patches++;
  }

  while (cur_rank < redist->size) {
    cur_rank++;
    recv_patches_by_rank[cur_rank] = &recv->recv_patches[recv->n_recv_patches];
    //mprintf("rank %d patches start at %d\n", cur_rank, recv->n_recv_patches);
  }

  int n_peers = 0;
  recv->buf_size = 0;
  for (int rank = 0; rank < redist->size; rank++) {
    struct mrc_redist_block *begin = recv_patches_by_rank[rank];
    struct mrc_redist_block *end   = recv_patches_by_rank[rank+1];

    // don't communicate with same proc, and skip procs that have no intersection
    if (rank == redist->rank || begin == end) {
      continue;
    }

    //mprintf("peer rank %d # = %ld\n", rank, end - begin);
    n_peers++;
  }
  //mprintf("n_peers %d\n", n_peers);

  recv->peers_begin = calloc(n_peers, sizeof(*recv->peers_begin));
  recv->peers_end = recv->peers_begin + n_peers;
  recv->reqs = calloc(n_peers, sizeof(*recv->reqs));

  struct mrc_redist_peer *peer = recv->peers_begin;
  for (int rank = 0; rank < redist->size; rank++) {
    struct mrc_redist_block *begin = recv_patches_by_rank[rank];
    struct mrc_redist_block *end   = recv_patches_by_rank[rank+1];

    // don't communicate with same proc, and skip procs that have no intersection
    if (rank == redist->rank || begin == end) {
      continue;
    }

    peer->rank = rank;
    peer->blocks_begin = begin;
    peer->blocks_end = end;

    // allocate buffer
    peer->buf_size = 0;
    for (struct mrc_redist_block *recv_patch = peer->blocks_begin; recv_patch < peer->blocks_end; recv_patch++) {
      int *ilo = recv_patch->ilo, *ihi = recv_patch->ihi;
      peer->buf_size += (size_t) (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);
    }
    
    recv->buf_size += peer->buf_size;
    recv->cnts[peer->rank] = peer->buf_size;
    peer++;
  }
  free(recv_patches_by_rank);

  // prefix sum to get displacements
  int last = 0;
  for (int r = 0; r < redist->size; r++) {
    recv->disps[r] = last;
    last += recv->cnts[r];
  }
  assert(last == recv->buf_size);
  
  // alloc aggregate recv buffers
  recv->buf = malloc(recv->buf_size * size_of_type);
  assert(recv->buf);
  int off = 0;
  for (struct mrc_redist_peer *peer = recv->peers_begin; peer != recv->peers_end; peer++) {
    peer->off = off;
    off += peer->buf_size;
  }

#if 0
  size_t g_data[2], data[2] = { recv->peers_end - recv->peers_begin, recv->buf_size };
  MPI_Reduce(data, g_data, 2, MPI_LONG, MPI_SUM, 0, redist->comm_writers);
  mpi_printf(redist->comm_writers, "avg recv buf size %ld (avg n_peers %ld)\n",
	     g_data[1] / redist->nr_writers, g_data[0] / redist->nr_writers);
#endif
}

// ----------------------------------------------------------------------
// mrc_redist_write_recv_post

static void
mrc_redist_write_recv_post(struct mrc_redist *redist, struct mrc_ndarray *nd,
			  struct mrc_fld *m3, int m)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;

  if (!redist->is_writer) {
    return;
  }
  
  for (struct mrc_redist_peer *peer = recv->peers_begin; peer < recv->peers_end; peer++) {
    // skip local patches
    if (peer->rank == redist->rank) {
      continue;
    }

    switch (mrc_fld_data_type(m3)) {
    case MRC_NT_FLOAT: {
      float *buf = recv->buf + peer->off * nd->size_of_type;
      for (struct mrc_redist_block *recv_patch = peer->blocks_begin; recv_patch < peer->blocks_end; recv_patch++) {
	int *ilo = recv_patch->ilo, *ihi = recv_patch->ihi;
	BUFLOOP(ix, iy, iz, ilo, ihi) {
	  MRC_S3(nd, ix,iy,iz) = *buf++;
	} BUFLOOP_END;
      }
      break;
    }
    case MRC_NT_DOUBLE: {
      double *buf = recv->buf + peer->off * nd->size_of_type;
      for (struct mrc_redist_block *recv_patch = peer->blocks_begin; recv_patch < peer->blocks_end; recv_patch++) {
	int *ilo = recv_patch->ilo, *ihi = recv_patch->ihi;
	BUFLOOP(ix, iy, iz, ilo, ihi) {
	  MRC_D3(nd, ix,iy,iz) = *buf++;
	} BUFLOOP_END;
      }
      break;
    }
    case MRC_NT_INT: {
      int *buf = recv->buf + peer->off * nd->size_of_type;
      for (struct mrc_redist_block *recv_patch = peer->blocks_begin; recv_patch < peer->blocks_end; recv_patch++) {
	int *ilo = recv_patch->ilo, *ihi = recv_patch->ihi;
	BUFLOOP(ix, iy, iz, ilo, ihi) {
	  MRC_I3(nd, ix,iy,iz) = *buf++;
	} BUFLOOP_END;
      }
      break;
    }
    default:
      assert(0);
    }    
  }
}

// ----------------------------------------------------------------------
// mrc_redist_write_destroy

static void
mrc_redist_write_destroy(struct mrc_redist *redist)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;

  free(recv->reqs);
  free(recv->disps);
  free(recv->buf);
  free(recv->peers_begin);

  free(recv->recv_patches);
}

// ----------------------------------------------------------------------
// mrc_redist_write_comm_local

static void
mrc_redist_write_comm_local(struct mrc_redist *redist, struct mrc_ndarray *nd,
			    struct mrc_fld *m3, int m)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->_domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch *patch = &patches[p];
    int *off = patch->off, *ldims = patch->ldims;

    int ilo[3], ihi[3];
    bool has_intersection =
      find_intersection(ilo, ihi, off, ldims, mrc_ndarray_offs(nd), mrc_ndarray_dims(nd));
    if (!has_intersection) {
      continue;
    }
    switch (mrc_fld_data_type(m3)) {
    case MRC_NT_FLOAT:
    {
      BUFLOOP(ix, iy, iz, ilo, ihi) {
        MRC_S3(nd, ix,iy,iz) =
          MRC_S5(m3, ix - off[0], iy - off[1], iz - off[2], m, p);
      } BUFLOOP_END
      break;
    }
    case MRC_NT_DOUBLE:
    {
      BUFLOOP(ix, iy, iz, ilo, ihi) {
        MRC_D3(nd, ix,iy,iz) =
          MRC_D5(m3, ix - off[0], iy - off[1], iz - off[2], m, p);
      } BUFLOOP_END
      break;     
    }
    case MRC_NT_INT:
    {
      BUFLOOP(ix, iy, iz, ilo, ihi) {
        MRC_I3(nd, ix,iy,iz) =
          MRC_I5(m3, ix - off[0], iy - off[1], iz - off[2], m, p);
      } BUFLOOP_END
      break;
    }
    default:
    {
      	assert(0);
    }
    }    
  }
}

// ----------------------------------------------------------------------
// mrc_redist_get_ndarray

struct mrc_ndarray *
mrc_redist_get_ndarray(struct mrc_redist *redist, struct mrc_fld *m3)
{
  mrc_redist_write_send_init(redist, m3->_nd->size_of_type);
  mrc_redist_write_recv_init(redist, m3->_nd->size_of_type);

  if (!redist->is_writer) {
    return NULL;
  }
  
  struct mrc_ndarray *nd = mrc_ndarray_create(redist->comm_writers);

  int writer;
  MPI_Comm_rank(redist->comm_writers, &writer);
  int writer_dims[3], writer_off[3];
  mrc_redist_writer_offs_dims(redist, writer, writer_off, writer_dims);
#if 0
  mprintf("writer_off %d %d %d dims %d %d %d\n",
	  writer_off[0], writer_off[1], writer_off[2],
	  writer_dims[0], writer_dims[1], writer_dims[2]);
#endif
  
  mrc_ndarray_set_param_int_array(nd, "dims", 3, writer_dims);
  mrc_ndarray_set_param_int_array(nd, "offs", 3, writer_off);
  
  switch (mrc_fld_data_type(m3)) {
  case MRC_NT_FLOAT: mrc_ndarray_set_type(nd, "float"); break;
  case MRC_NT_DOUBLE: mrc_ndarray_set_type(nd, "double"); break;
  case MRC_NT_INT: mrc_ndarray_set_type(nd, "int"); break;
  default: assert(0);
  }
  mrc_ndarray_setup(nd);

  mrc_redist_write_recv_init2(redist, nd, m3->_nd->size_of_type);

  return nd;
}

// ----------------------------------------------------------------------
// mrc_redist_put_ndarray

void
mrc_redist_put_ndarray(struct mrc_redist *redist, struct mrc_ndarray *nd)
{
  mrc_redist_write_send_destroy(redist);

  if (redist->is_writer) {
    mrc_redist_write_destroy(redist);
    mrc_ndarray_destroy(nd);
  }
}

// ----------------------------------------------------------------------
// mrc_redist_write_begin

static void
mrc_redist_write_begin(struct mrc_redist *redist, struct mrc_ndarray *nd,
		       struct mrc_fld *m3, int m)
{
  mrc_redist_write_send_prep(redist, m3, m);

  struct mrc_redist_write_recv *recv = &redist->write_recv;
  struct mrc_redist_write_send *send = &redist->write_send;
  MPI_Datatype mpi_dtype = to_mpi_datatype(mrc_fld_data_type(m3));
  
  if (redist->is_writer) {
    for (struct mrc_redist_peer *peer = recv->peers_begin; peer < recv->peers_end; peer++) {
      //mprintf("recv_begin: Irecv cnt %ld from %d\n", peer->buf_size, peer->rank);
      assert(peer->off == recv->disps[peer->rank]);
      MPI_Irecv(recv->buf + recv->disps[peer->rank] * nd->size_of_type, peer->buf_size,
		mpi_dtype, peer->rank, 0x1000, redist->comm,
		&recv->reqs[peer - recv->peers_begin]);
    }
  }

  for (struct mrc_redist_peer *w = send->peers_begin; w != send->peers_end; w++) {
    //mprintf("send_begin: Isend cnt %ld to %d\n", w->buf_size, w->rank);
    MPI_Isend(send->buf + send->disps[w->rank] * m3->_nd->size_of_type, w->buf_size,
	      mpi_dtype, w->rank, 0x1000, redist->comm,
	      &send->reqs[w - send->peers_begin]);
  }
}

// ----------------------------------------------------------------------
// mrc_redist_write_end

static void
mrc_redist_write_end(struct mrc_redist *redist, struct mrc_ndarray *nd,
		     struct mrc_fld *m3, int m)
{
  struct mrc_redist_write_recv *recv = &redist->write_recv;
  struct mrc_redist_write_send *send = &redist->write_send;

  if (redist->is_writer) {
    //mprintf("recv_end: Waitall cnt = %d\n", recv->peers_end - recv->peers_begin);
    MPI_Waitall(recv->peers_end - recv->peers_begin, recv->reqs, MPI_STATUSES_IGNORE);
  }

  int n_peers = send->peers_end - send->peers_begin;
  //mprintf("send_end: Waitall cnt = %d\n", n_peers);
  MPI_Waitall(n_peers, send->reqs, MPI_STATUSES_IGNORE);

  mrc_redist_write_recv_post(redist, nd, m3, m);
}

// ----------------------------------------------------------------------
// mrc_redist_write_begin2

static void
mrc_redist_write_begin2(struct mrc_redist *redist, struct mrc_ndarray *nd,
			struct mrc_fld *m3, int m)
{
  mrc_redist_write_send_prep(redist, m3, m);

  struct mrc_redist_write_recv *recv = &redist->write_recv;
  struct mrc_redist_write_send *send = &redist->write_send;
  MPI_Datatype mpi_dtype = to_mpi_datatype(mrc_fld_data_type(m3));

#if 0
  for (int r = 0; r < redist->size; r++) {
    MPI_Barrier(redist->comm);
    if (r == redist->rank) {
      mprintf("send disp/cnts:");
      for (int r = 0; r < redist->size; r++) {
	printf(" (%d) %d/%d", r, send->disps[r], send->cnts[r]);
      }
      printf("\n");
      mprintf("recv disp/cnts:");
      for (int r = 0; r < redist->size; r++) {
	printf(" (%d) %d/%d", r, recv->disps[r], recv->cnts[r]);
      }
      printf("\n");
    }
  }
#endif
  
  MPI_Barrier(redist->comm);
  MPI_Alltoallv(send->buf, send->cnts, send->disps, mpi_dtype,
  		recv->buf, recv->cnts, recv->disps, mpi_dtype, redist->comm);
}

// ----------------------------------------------------------------------
// mrc_redist_write_end2

static void
mrc_redist_write_end2(struct mrc_redist *redist, struct mrc_ndarray *nd,
		      struct mrc_fld *m3, int m)
{
  mrc_redist_write_recv_post(redist, nd, m3, m);
}

// ----------------------------------------------------------------------

void
mrc_redist_run(struct mrc_redist *redist, struct mrc_ndarray *nd,
	       struct mrc_fld *m3, int m)
{
  mrc_redist_write_begin2(redist, nd, m3, m);

  if (redist->is_writer) {
    mrc_redist_write_comm_local(redist, nd, m3, m);
  }

  mrc_redist_write_end2(redist, nd, m3, m);
}



