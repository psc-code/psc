
#ifndef DDC_PARTICLES_HXX
#define DDC_PARTICLES_HXX

#include <mrc_ddc.h>
#include <mrc_domain.h>

// ======================================================================
// ddc_particles

#define N_DIR (27)

template<typename MP>
struct ddc_particles
{
  using Mparticles = MP;
  using particle_t = typename Mparticles::particle_t;
  using buf_t = typename Mparticles::buf_t;
  using real_t = typename Mparticles::real_t;
  
  ddc_particles(const Grid_t& grid);

  void comm(std::vector<buf_t*>& bufs);

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
    std::vector<dsend_entry> send_entry;
    std::vector<int> send_cnts;
    int n_send_entries;
    int n_send;
    
    std::vector<drecv_entry> recv_entry;
    std::vector<int> recv_cnts;
    int n_recv_entries;
    int n_recv;
    
    int rank;
  };

  struct dnei {
    buf_t send_buf;
    int n_recv;
    int rank;
    int patch;
  };

  struct patch {
    dnei nei[N_DIR];
    int n_recv;
  };
  
  int nr_patches;
  std::vector<patch> patches_;
  std::vector<ddcp_info_by_rank> cinfo_; // compressed info
  int n_ranks;
  std::vector<MPI_Request> send_reqs_;
  std::vector<MPI_Request> recv_reqs_;
};

// ----------------------------------------------------------------------
// ctor

template<typename MP>
inline ddc_particles<MP>::ddc_particles(const Grid_t& grid)
{
  std::memset(this, 0, sizeof(*this));

  nr_patches = grid.n_patches();
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
	    grid.neighborRankPatch(p, dir, &nei->rank, &nei->patch);
	  }
	}
      }
    }
  }

  MPI_Comm comm = grid.comm();
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  std::vector<ddcp_info_by_rank> info(size);

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
	  info[nei->rank].n_recv_entries++;
	  info[nei->rank].rank = nei->rank;
	}
      }
    }
  }

  // alloc recv_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      info[r].recv_entry.resize(info[r].n_recv_entries);
      info[r].recv_cnts.resize(info[r].n_recv_entries);
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
	  info[nei->rank].n_send_entries++;
	}
      }
    }
  }

  // alloc send_entries
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      info[r].send_entry.reserve(info[r].n_send_entries);
      info[r].send_cnts.resize(info[r].n_send_entries);
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
	  auto se = dsend_entry{};
	  se.patch = p;
	  se.nei_patch = nei->patch;
	  se.dir1 = dir1;
	  se.dir1neg = dir1neg;
	  info[nei->rank].send_entry.push_back(se);
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
  n_ranks = n_send_ranks;

  send_reqs_.resize(n_ranks);
  recv_reqs_.resize(n_ranks);

  n_recv_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_recv_entries) {
      MPI_Irecv(info[r].recv_entry.data(),
		sizeof(drecv_entry) / sizeof(int) * info[r].n_recv_entries,
		MPI_INT, r, 111, comm, &recv_reqs_[n_recv_ranks++]);
    }
  }  

  n_send_ranks = 0;
  for (int r = 0; r < size; r++) {
    if (info[r].n_send_entries) {
      MPI_Isend(info[r].send_entry.data(),
		sizeof(dsend_entry) / sizeof(int) * info[r].n_send_entries,
		MPI_INT, r, 111, comm, &send_reqs_[n_send_ranks++]);
    }
  }  

  cinfo_.reserve(n_ranks);
  for (auto& item : info) {
    if (item.n_recv_entries) {
      assert(item.n_send_entries);
      cinfo_.push_back(std::move(item));
    }
  }
  assert(cinfo_.size() == n_ranks);
}

// ----------------------------------------------------------------------
// comm
//
// OPT: could use MPI_Waitany?
// OPT: overall more async
// OPT: 1d instead of 3d loops
// OPT: make the status buffers only as large as needed?

template<typename MP>
inline void ddc_particles<MP>::comm(std::vector<buf_t*>& bufs)
{
  using iterator_t = typename buf_t::iterator;
  
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // FIXME, this is assuming our struct is equiv to an array of real_type
  assert(sizeof(particle_t) % sizeof(real_t) == 0);
  int sz = sizeof(particle_t) / sizeof(real_t);
  int dir[3];

  for (int r = 0; r < n_ranks; r++) {
    MPI_Irecv(cinfo_[r].recv_cnts.data(), cinfo_[r].n_recv_entries,
	      MPI_INT, cinfo_[r].rank, 222, comm, &recv_reqs_[r]);
  }  

  for (int r = 0; r < n_ranks; r++) {
    cinfo_[r].n_send = 0;
    for (int i = 0; i < cinfo_[r].n_send_entries; i++) {
      dsend_entry *se = &cinfo_[r].send_entry[i];
      patch *patch = &patches_[se->patch];
      dnei *nei = &patch->nei[se->dir1];
      unsigned int n_send = nei->send_buf.size();
      cinfo_[r].send_cnts[i] = n_send;
      cinfo_[r].n_send += n_send;
    }
    MPI_Isend(cinfo_[r].send_cnts.data(), cinfo_[r].n_send_entries,
	      MPI_INT, cinfo_[r].rank, 222, comm, &send_reqs_[r]);
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

  MPI_Waitall(n_ranks, recv_reqs_.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(n_ranks, send_reqs_.data(), MPI_STATUSES_IGNORE);

  MPI_Datatype mpi_dtype = Mparticles_traits<Mparticles>::mpi_dtype();

  // add remote # particles
  int n_send = 0, n_recv = 0;

  for (int r = 0; r < n_ranks; r++) {
    cinfo_[r].n_recv = 0;
    for (int i = 0; i < cinfo_[r].n_recv_entries; i++) {
      drecv_entry *re = &cinfo_[r].recv_entry[i];
      patch *patch = &patches_[re->patch];
      patch->n_recv += cinfo_[r].recv_cnts[i];
      cinfo_[r].n_recv += cinfo_[r].recv_cnts[i];
    }
    n_send += cinfo_[r].n_send;
    n_recv += cinfo_[r].n_recv;
  }

  // post sends
  buf_t send_buf;
  send_buf.resize(n_send);
  iterator_t it = send_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    if (cinfo_[r].n_send == 0)
      continue;

    iterator_t it0 = it;
    for (int i = 0; i < cinfo_[r].n_send_entries; i++) {
      dsend_entry *se = &cinfo_[r].send_entry[i];
      patch *patch = &patches_[se->patch];
      buf_t *send_buf_nei = &patch->nei[se->dir1].send_buf;
      std::copy(send_buf_nei->begin(), send_buf_nei->end(), it);
      it += send_buf_nei->size();
    }
    MPI_Isend(&*it0, sz * cinfo_[r].n_send, mpi_dtype,
	      cinfo_[r].rank, 1, comm, &send_reqs_[r]);
  }
  assert(it == send_buf.begin() + n_send);

  // post receives
  buf_t recv_buf;
  recv_buf.resize(n_recv);
  it = recv_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    if (cinfo_[r].n_recv == 0)
      continue;

    MPI_Irecv(&*it, sz * cinfo_[r].n_recv, mpi_dtype,
	      cinfo_[r].rank, 1, comm, &recv_reqs_[r]);
    it += cinfo_[r].n_recv;
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
    auto& buf = *bufs[p];
    int size = buf.size();
    buf.reserve(size + patch->n_recv);
    // this is dangerous: we keep using the iterator, knowing that
    // it won't become invalid due to a realloc since we reserved enough space...
    it_recv[p] = buf.end();
    buf.resize(size + patch->n_recv);
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
	  buf_t *nei_send_buf = &patches_[nei->patch].nei[dir1neg].send_buf;

	  std::copy(nei_send_buf->begin(), nei_send_buf->end(), it_recv[p]);
	  it_recv[p] += nei_send_buf->size();
	}
      }
    }
  }

  MPI_Waitall(n_ranks, recv_reqs_.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(n_ranks, send_reqs_.data(), MPI_STATUSES_IGNORE);

  // copy received particles into right place

  it = recv_buf.begin();
  for (int r = 0; r < n_ranks; r++) {
    for (int i = 0; i < cinfo_[r].n_recv_entries; i++) {
      drecv_entry *re = &cinfo_[r].recv_entry[i];
      std::copy(it, it + cinfo_[r].recv_cnts[i], it_recv[re->patch]);
      it += cinfo_[r].recv_cnts[i];
      it_recv[re->patch] += cinfo_[r].recv_cnts[i];
    }
  }
  assert(it == recv_buf.begin() + n_recv);

  for (int p = 0; p < nr_patches; p++) {
    assert(it_recv[p] == bufs[p]->end());
  }
  
  delete[] it_recv;
}

#endif

