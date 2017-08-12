
#define DDCP_TYPE_COMMON     1
#define DDCP_TYPE_COMMON2    2
#define DDCP_TYPE_COMMON_OMP 3
#define DDCP_TYPE_CUDA       4

#include <mrc_ddc.h>

#include <string.h>

// ----------------------------------------------------------------------
// ddcp_buf_t

#if DDCP_TYPE == DDCP_TYPE_COMMON || DDCP_TYPE == DDCP_TYPE_COMMON2 || DDCP_TYPE == DDCP_TYPE_COMMON_OMP

static void
ddcp_particles_realloc(struct psc_mparticles *mprts, int p, int new_n_particles)
{
  mparticles_patch_reserve(mprts, p, new_n_particles);
}

static particle_t *
ddcp_particles_get_addr(struct psc_mparticles *mprts, int p, int n)
{
  particle_range_t prts = particle_range_mprts(mprts, p);
  return particle_iter_at(prts.begin, n);
}

#elif DDCP_TYPE == DDCP_TYPE_CUDA

static void
ddcp_particles_realloc(struct psc_mparticles *mprts, int p, int new_n_particles)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cmprts->bnd.bpatch[p].prts = realloc(cmprts->bnd.bpatch[p].prts, new_n_particles * sizeof(*cmprts->bnd.bpatch[p].prts));
}

static particle_t *
ddcp_particles_get_addr(struct psc_mparticles *mprts, int p, int n)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  return &cmprts->bnd.bpatch[p].prts[n];
}

#endif

static void ddcp_particles_realloc(struct psc_mparticles *mprts, int p, int new_n_particles);
static particle_t *ddcp_particles_get_addr(struct psc_mparticles *mprts, int p, int n);

#define N_DIR (27)

typedef struct {
  struct psc_mparticles *m_mprts;
  int m_p;
} ddcp_buf_t;

static void
ddcp_buf_ctor(ddcp_buf_t *buf, struct psc_mparticles *mprts, int p)
{
  buf->m_mprts = mprts;
  buf->m_p = p;
}

static void
ddcp_buf_dtor(ddcp_buf_t *buf)
{
}

static particle_t *
ddcp_buf_at(ddcp_buf_t *buf, int n)
{
  return ddcp_particles_get_addr(buf->m_mprts, buf->m_p, n);
}

static void
ddcp_buf_reserve(ddcp_buf_t *buf, int n)
{
  ddcp_particles_realloc(buf->m_mprts, buf->m_p, n);
}

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
  int head;
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

  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp_buf_ctor(&ddcp->patches[p].buf, mprts, p);
  }
  
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
    ddcp_buf_reserve(&patch->buf, patch->head + patch->n_recv);
  }

  // leave room for receives (FIXME? just change order)
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      patch->head += cinfo[r].recv_cnts[i];
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
	  particle_t *patch_it = ddcp_buf_at(&patch->buf, patch->head);
	  particle_buf_t *nei_send_buf = &ddcp->patches[nei->patch].nei[dir1neg].send_buf;
	  particle_buf_copy(particle_buf_begin(nei_send_buf), particle_buf_end(nei_send_buf),
			    patch_it);
	  patch->head += particle_buf_size(nei_send_buf);
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
    old_head[p] = patch->head;
    patch->head -= patch->n_recv;
  }

  it = particle_buf_begin(&recv_buf);
  for (int r = 0; r < ddcp->n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      struct ddcp_recv_entry *re = &cinfo[r].recv_entry[i];
      struct ddcp_patch *patch = &ddcp->patches[re->patch];
      particle_t *patch_it = ddcp_buf_at(&patch->buf, patch->head);
      particle_buf_copy(it, it + cinfo[r].recv_cnts[i], patch_it);
      patch->head += cinfo[r].recv_cnts[i];
      it += cinfo[r].recv_cnts[i];
    }
  }
  assert(it == particle_buf_begin(&recv_buf) + n_recv);

  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    patch->head = old_head[p];
  }
  free(old_head);

  particle_buf_dtor(&send_buf);
  particle_buf_dtor(&recv_buf);

  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp_buf_dtor(&ddcp->patches[p].buf);
  }
}

