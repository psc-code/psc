
#ifndef DDC_PARTICLES_HXX
#define DDC_PARTICLES_HXX

#include <mrc_ddc.h>

// ======================================================================
// ddc_particles

#define N_DIR (27)

template<typename MP>
struct ddc_particles
{
  using Mparticles = MP;
  using particle_t = typename Mparticles::particle_t;
  using particle_buf_t = typename Mparticles::buf_t;
  using real_t = typename Mparticles::real_t;
  
  ddc_particles(struct mrc_domain *domain);
  ~ddc_particles();

  void comm();

  struct dsend_entry {
    int patch; // source patch (source rank is this rank)
    int nei_patch; // target patch (target rank is index in send_entry)
    int dir1;  // direction
    int dir1neg;
  };
  
  struct drecv_entry { // needs to be same as send_entry with different order!
    int nei_patch;
    int patch;
    int dir1neg;
    int dir1;
  };
  
  struct ddcp_info_by_rank {
    dsend_entry *send_entry;
    int *send_cnts;
    int n_send_entries;
    int n_send;
    
    drecv_entry *recv_entry;
    int *recv_cnts;
    int n_recv_entries;
    int n_recv;
    
    int rank;
  };

  struct dnei {
    particle_buf_t send_buf;
    int n_recv;
    int rank;
    int patch;
  };

  struct patch {
    particle_buf_t *m_buf;
    dnei nei[N_DIR];
    int n_recv;
  };
  
  int nr_patches;
  std::vector<patch> patches_;
  std::vector<ddcp_info_by_rank> by_rank_;
  struct ddcp_info_by_rank *cinfo; // compressed info
  int n_ranks;
  MPI_Request *send_reqs;
  MPI_Request *recv_reqs;

  struct mrc_domain *domain;
};

// ----------------------------------------------------------------------
// ctor

template<typename MP>
inline ddc_particles<MP>::ddc_particles(struct mrc_domain *_domain)
{
  std::memset(this, 0, sizeof(*this));

  domain = _domain;
  mrc_domain_get_patches(domain, &nr_patches);
  patches_.resize(nr_patches);
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  dnei *nei = &patch->nei[dir1];

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

  by_rank_.resize(size);

  int dir[3];

  // count how many recv_entries per rank
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  dnei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  by_rank_[nei->rank].n_recv_entries++;
	}
      }
    }
  }

  // alloc recv_entries
  for (int r = 0; r < size; r++) {
    if (by_rank_[r].n_recv_entries) {
      by_rank_[r].recv_entry =
	(drecv_entry *) malloc(by_rank_[r].n_recv_entries * sizeof(*by_rank_[r].recv_entry));
      by_rank_[r].recv_cnts =
	(int *) malloc(by_rank_[r].n_recv_entries * sizeof(*by_rank_[r].recv_cnts));
    }
  }

#if 0
  // set up recv_entries
  for (int p = 0; p < nr_patches; p++) {
    ddc_particles::patch *patch = &patches[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  dnei *nei = &patch->nei[dir1];
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
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  dnei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  by_rank_[nei->rank].n_send_entries++;
	}
      }
    }
  }

  // alloc send_entries
  for (int r = 0; r < size; r++) {
    if (by_rank_[r].n_send_entries) {
      by_rank_[r].send_entry = (dsend_entry *)
	malloc(by_rank_[r].n_send_entries * sizeof(*by_rank_[r].send_entry));
      by_rank_[r].send_cnts = (int *)
	malloc(by_rank_[r].n_send_entries * sizeof(*by_rank_[r].send_cnts));
      by_rank_[r].n_send_entries = 0;
    }
  }

  // set up send_entries
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dirneg[3] = { -dir[0], -dir[1], -dir[2] };
	  int dir1neg = mrc_ddc_dir2idx(dirneg);
	  dnei *nei = &patch->nei[dir1];
	  if (nei->rank < 0 || nei->rank == rank) {
	    continue;
	  }
	  dsend_entry *se =
	    &by_rank_[nei->rank].send_entry[by_rank_[nei->rank].n_send_entries++];
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
    if (by_rank_[r].n_recv_entries) {
      n_recv_ranks++;
    }
    if (by_rank_[r].n_send_entries) {
      n_send_ranks++;
    }
  }

  assert(n_send_ranks == n_recv_ranks);
  n_ranks = n_send_ranks;

  send_reqs = (MPI_Request *) malloc(n_ranks * sizeof(*send_reqs));
  recv_reqs = (MPI_Request *) malloc(n_ranks * sizeof(*recv_reqs));

  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (by_rank_[r].n_recv_entries) {
      MPI_Irecv(by_rank_[r].recv_entry,
		sizeof(drecv_entry) / sizeof(int) * by_rank_[r].n_recv_entries,
		MPI_INT, r, 111, comm, &recv_reqs[n_recv_ranks++]);
    }
  }  

  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (by_rank_[r].n_send_entries) {
      MPI_Isend(by_rank_[r].send_entry,
		sizeof(dsend_entry) / sizeof(int) * by_rank_[r].n_send_entries,
		MPI_INT, r, 111, comm, &send_reqs[n_send_ranks++]);
    }
  }  

  // FIXME / OPT, we're copying alloc'd pointers over,
  // fragile, though correct. info could be free'd here or even
  // earlier
  cinfo = (struct ddcp_info_by_rank *) malloc(n_ranks * sizeof(*cinfo));
  int i = 0;
  for (int r = 0; r < size; r++) {
    if (by_rank_[r].n_recv_entries) {
      assert(by_rank_[r].n_send_entries);
      cinfo[i] = by_rank_[r];
      cinfo[i].rank = r;
      i++;
    }
  }
  assert(i == n_ranks);
}

// ----------------------------------------------------------------------
// dtor

template<typename MP>
inline ddc_particles<MP>::~ddc_particles()
{
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    int dir[3];
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	}
      }
    }
  }

  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int size;
  MPI_Comm_size(comm, &size);

  for (int r = 0; r < size; r++) {
    free(by_rank_[r].send_entry);
    free(by_rank_[r].recv_entry);
    free(by_rank_[r].send_cnts);
    free(by_rank_[r].recv_cnts);
  }
  free(cinfo);
  free(recv_reqs);
  free(send_reqs);
}

// ----------------------------------------------------------------------
// comm
//
// OPT: could use MPI_Waitany?
// OPT: overall more async
// OPT: 1d instead of 3d loops
// OPT: make the status buffers only as large as needed?

template<typename MP>
inline void ddc_particles<MP>::comm()
{
  using iterator_t = typename particle_buf_t::iterator;
  
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // FIXME, this is assuming our struct is equiv to an array of real_type
  assert(sizeof(particle_t) % sizeof(real_t) == 0);
  int sz = sizeof(particle_t) / sizeof(real_t);
  int dir[3];

  for (int r = 0; r < n_ranks; r++) {
    MPI_Irecv(cinfo[r].recv_cnts, cinfo[r].n_recv_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &recv_reqs[r]);
  }  

  for (int r = 0; r < n_ranks; r++) {
    cinfo[r].n_send = 0;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      dsend_entry *se = &cinfo[r].send_entry[i];
      patch *patch = &patches_[se->patch];
      dnei *nei = &patch->nei[se->dir1];
      unsigned int n_send = nei->send_buf.size();
      cinfo[r].send_cnts[i] = n_send;
      cinfo[r].n_send += n_send;
    }
    MPI_Isend(cinfo[r].send_cnts, cinfo[r].n_send_entries,
	      MPI_INT, cinfo[r].rank, 222, comm, &send_reqs[r]);
  }  

  // overlap: count local # particles
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];
    patch->n_recv = 0;
    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  dnei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  int dirneg[3] = { -dir[0], -dir[1], -dir[2] };
	  int dir1neg = mrc_ddc_dir2idx(dirneg);
	  dnei *nei_send = &patches_[nei->patch].nei[dir1neg];
	  patch->n_recv += nei_send->send_buf.size();
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

  MPI_Waitall(n_ranks, recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(n_ranks, send_reqs, MPI_STATUSES_IGNORE);

  MPI_Datatype mpi_dtype = Mparticles_traits<Mparticles>::mpi_dtype();

  // add remote # particles
  int n_send = 0, n_recv = 0;

  for (int r = 0; r < n_ranks; r++) {
    cinfo[r].n_recv = 0;
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      drecv_entry *re = &cinfo[r].recv_entry[i];
      patch *patch = &patches_[re->patch];
      patch->n_recv += cinfo[r].recv_cnts[i];
      cinfo[r].n_recv += cinfo[r].recv_cnts[i];
    }
    n_send += cinfo[r].n_send;
    n_recv += cinfo[r].n_recv;
  }

  // post sends
  particle_buf_t send_buf;
  send_buf.resize(n_send);
  iterator_t it = send_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    if (cinfo[r].n_send == 0)
      continue;

    iterator_t it0 = it;
    for (int i = 0; i < cinfo[r].n_send_entries; i++) {
      dsend_entry *se = &cinfo[r].send_entry[i];
      patch *patch = &patches_[se->patch];
      particle_buf_t *send_buf_nei = &patch->nei[se->dir1].send_buf;
      std::copy(send_buf_nei->begin(), send_buf_nei->end(), it);
      it += send_buf_nei->size();
    }
    MPI_Isend(&*it0, sz * cinfo[r].n_send, mpi_dtype,
	      cinfo[r].rank, 1, comm, &send_reqs[r]);
  }
  assert(it == send_buf.begin() + n_send);

  // post receives
  particle_buf_t recv_buf;
  recv_buf.resize(n_recv);
  it = recv_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    if (cinfo[r].n_recv == 0)
      continue;

    MPI_Irecv(&*it, sz * cinfo[r].n_recv, mpi_dtype,
	      cinfo[r].rank, 1, comm, &recv_reqs[r]);
    it += cinfo[r].n_recv;
  }
  assert(it == recv_buf.begin() + n_recv);

  // leave room for receives (FIXME? just change order)
  // each patch's array looks like:
  // [........|.........|...]
  // ---------                    particles remaining in this patch
  //          --------------      new particles go here (# = patch->n_recvs)
  //          ----------          locally exchanged particles go here
  //                    ----      remote particles go here
  iterator_t *it_recv = new iterator_t[nr_patches];

  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];
    int size = patch->m_buf->size();
    patch->m_buf->reserve(size + patch->n_recv);
    // this is dangerous: we keep using the iterator, knowing that
    // it won't become invalid due to a realloc since we reserved enough space...
    it_recv[p] = patch->m_buf->end();
    patch->m_buf->resize(size + patch->n_recv);
  }

  // overlap: copy particles from local proc to the end of recv range
  for (int p = 0; p < nr_patches; p++) {
    patch *patch = &patches_[p];

    for (dir[2] = -1; dir[2] <= 1; dir[2]++) {
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  int dir1 = mrc_ddc_dir2idx(dir);
	  int dirneg[3] = { -dir[0], -dir[1], -dir[2] };
	  int dir1neg = mrc_ddc_dir2idx(dirneg);
	  dnei *nei = &patch->nei[dir1];
	  if (nei->rank != rank) {
	    continue;
	  }
	  particle_buf_t *nei_send_buf = &patches_[nei->patch].nei[dir1neg].send_buf;

	  std::copy(nei_send_buf->begin(), nei_send_buf->end(), it_recv[p]);
	  it_recv[p] += nei_send_buf->size();
	}
      }
    }
  }

  MPI_Waitall(n_ranks, recv_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(n_ranks, send_reqs, MPI_STATUSES_IGNORE);

  // copy received particles into right place

  it = recv_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    for (int i = 0; i < cinfo[r].n_recv_entries; i++) {
      drecv_entry *re = &cinfo[r].recv_entry[i];
      std::copy(it, it + cinfo[r].recv_cnts[i], it_recv[re->patch]);
      it += cinfo[r].recv_cnts[i];
      it_recv[re->patch] += cinfo[r].recv_cnts[i];
    }
  }
  assert(it == recv_buf.begin() + n_recv);

  for (int p = 0; p < nr_patches; p++) {
    assert(it_recv[p] == patches_[p].m_buf->end());
  }
  
  delete[] it_recv;
}

#endif

