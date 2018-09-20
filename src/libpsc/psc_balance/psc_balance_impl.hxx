
#pragma once

#include "fields.hxx"
#include "balance.hxx"
#include "bnd_particles.hxx"
#include "bnd.hxx"

#include <mrc_profile.h>
#include <string.h>

extern double *psc_balance_comp_time_by_patch;

static double
capability_default(int p)
{
  return 1.;
}

static double _mrc_unused
capability_jaguar(int p)
{
  if (p % 16 == 0) {
    return 16.;
  } else {
    return 1.;
  }
}

// ======================================================================
// Communicate

struct send_info {
  int rank;
  int patch;
};

struct recv_info {
  int rank;
  int patch;
};

struct by_ri {
  int rank;
  int nr_patches;
  int *pi_to_patch;
};

struct communicate_ctx {
  MPI_Comm comm;
  int mpi_rank;
  int mpi_size;
  int nr_patches_old;
  int nr_patches_new;
  struct send_info *send_info; // by old patch on this proc
  struct recv_info *recv_info; // by new patch on this proc

  int *send_rank_to_ri; // map from send target rank to contiguous "rank index"
  int *recv_rank_to_ri; // map from recv source rank to contiguous "rank index"

  int nr_send_ranks;
  struct by_ri *send_by_ri;

  int nr_recv_ranks;
  struct by_ri *recv_by_ri;

  communicate_ctx(const MrcDomain& domain_old, const MrcDomain& domain_new)
  {
    comm = domain_old.comm();
    MPI_Comm_rank(comm, &mpi_rank);
    MPI_Comm_size(comm, &mpi_size);

    nr_patches_old = domain_old.nPatches();
    nr_patches_new = domain_new.nPatches();

    send_info = (struct send_info *) calloc(nr_patches_old, sizeof(*send_info));
    recv_info = (struct recv_info *) calloc(nr_patches_new, sizeof(*recv_info));

    for (int p = 0; p < nr_patches_old; p++) {
      auto info_old = domain_old.localPatchInfo(p);
      auto info_new = domain_new.levelIdx3PatchInfo(info_old.level, info_old.idx3);
      send_info[p].rank  = info_new.rank;
      send_info[p].patch = info_new.patch;
    }
    for (int p = 0; p < nr_patches_new; p++) {
      auto info_new = domain_new.localPatchInfo(p);
      auto info_old = domain_old.levelIdx3PatchInfo(info_new.level, info_new.idx3);
      recv_info[p].rank  = info_old.rank;
      recv_info[p].patch = info_old.patch;
    }

    // maps rank <-> rank index

    send_rank_to_ri = (int *) malloc(mpi_size * sizeof(*send_rank_to_ri));
    recv_rank_to_ri = (int *) malloc(mpi_size * sizeof(*recv_rank_to_ri));
    for (int r = 0; r < mpi_size; r++) {
      send_rank_to_ri[r] = -1;
      recv_rank_to_ri[r] = -1;
    }

    nr_send_ranks = 0;
    nr_recv_ranks = 0;
    for (int p = 0; p < nr_patches_old; p++) {
      int send_rank = send_info[p].rank;
      if (send_rank >= 0) {
	if (send_rank_to_ri[send_rank] < 0) {
	  send_rank_to_ri[send_rank] = nr_send_ranks++;
	}
      }
    }
    for (int p = 0; p < nr_patches_new; p++) {
      int recv_rank = recv_info[p].rank;
      if (recv_rank >= 0) {
	if (recv_rank_to_ri[recv_rank] < 0) {
	  recv_rank_to_ri[recv_rank] = nr_recv_ranks++;
	}
      }
    }

    send_by_ri = (struct by_ri *) calloc(nr_send_ranks, sizeof(*send_by_ri));
    recv_by_ri = (struct by_ri *) calloc(nr_recv_ranks, sizeof(*recv_by_ri));

    for (int p = 0; p < nr_patches_old; p++) {
      int send_rank = send_info[p].rank;
      if (send_rank >= 0) {
	send_by_ri[send_rank_to_ri[send_rank]].rank = send_rank;
      }
    }
    for (int p = 0; p < nr_patches_new; p++) {
      int recv_rank = recv_info[p].rank;
      if (recv_rank >= 0) {
	recv_by_ri[recv_rank_to_ri[recv_rank]].rank = recv_rank;
      }
    }

    /* for (int ri = 0; ri < nr_send_ranks; ri++) { */
    /*   mprintf("send -> %d (%d)\n", send_by_ri[ri].rank, ri); */
    /* } */

    // count number of patches sent to each rank

    for (int p = 0; p < nr_patches_old; p++) {
      int send_rank = send_info[p].rank;
      if (send_rank >= 0) {
	send_by_ri[send_rank_to_ri[send_rank]].nr_patches++;
      }
    }

    // map send patch index by ri back to local patch number

    for (int ri = 0; ri < nr_send_ranks; ri++) {
      send_by_ri[ri].pi_to_patch = (int *) calloc(send_by_ri[ri].nr_patches, sizeof(*send_by_ri[ri].pi_to_patch));
      send_by_ri[ri].nr_patches = 0;
    }

    for (int p = 0; p < nr_patches_old; p++) {
      int send_rank = send_info[p].rank;
      if (send_rank < 0) {
	continue;
      }
      int ri = send_rank_to_ri[send_rank];
      int pi = send_by_ri[ri].nr_patches++;
      send_by_ri[ri].pi_to_patch[pi] = p;
    }

    // count number of patches received from each rank

    for (int p = 0; p < nr_patches_new; p++) {
      int recv_rank = recv_info[p].rank;
      if (recv_rank >= 0) {
	int ri = recv_rank_to_ri[recv_rank];
	recv_by_ri[ri].nr_patches++;
      }
    }

    // map received patch index by ri back to local patch number

    for (int ri = 0; ri < nr_recv_ranks; ri++) {
      recv_by_ri[ri].pi_to_patch = (int *) calloc(recv_by_ri[ri].nr_patches, sizeof(*recv_by_ri[ri].pi_to_patch));
      recv_by_ri[ri].nr_patches = 0;
    }

    for (int p = 0; p < nr_patches_new; p++) {
      int recv_rank = recv_info[p].rank;
      if (recv_rank < 0) {
	continue;
      }
      int ri = recv_rank_to_ri[recv_rank];
      int pi = recv_by_ri[ri].nr_patches++;
      recv_by_ri[ri].pi_to_patch[pi] = p;
    }

  }

  ~communicate_ctx()
  {
    free(send_info);
    free(recv_info);
    
    free(send_rank_to_ri);
    free(recv_rank_to_ri);
    
    for (int ri = 0; ri < nr_send_ranks; ri++) {
      free(send_by_ri[ri].pi_to_patch);
    }
    free(send_by_ri);
    
    for (int ri = 0; ri < nr_recv_ranks; ri++) {
      free(recv_by_ri[ri].pi_to_patch);
    }
    free(recv_by_ri);
  }

  std::vector<uint> new_n_prts(std::vector<uint> n_prts_by_patch_old)
  {
    static int pr;
    if (!pr) {
      pr   = prof_register("comm nr prts", 1., 0, 0);
    }

    prof_start(pr);

    std::vector<uint> n_prts_by_patch_new(nr_patches_new);
    // post receives 

    MPI_Request *recv_reqs = (MPI_Request *) calloc(nr_patches_new, sizeof(*recv_reqs));
    int nr_recv_reqs = 0;

    int **nr_particles_recv_by_ri = (int **) calloc(nr_recv_ranks, sizeof(*nr_particles_recv_by_ri));
    for (int ri = 0; ri < nr_recv_ranks; ri++) {
      struct by_ri *recv = &recv_by_ri[ri];
      nr_particles_recv_by_ri[ri] = (int *) calloc(recv->nr_patches, sizeof(*nr_particles_recv_by_ri[ri]));
    
      if (recv->rank != mpi_rank) {
	//mprintf("recv <- %d (len %d)\n", r, nr_patches_recv_by_ri[ri]);
	MPI_Irecv(nr_particles_recv_by_ri[ri], recv->nr_patches, MPI_INT,
		  recv->rank, 10, comm, &recv_reqs[nr_recv_reqs++]);
      }
    }

    // post sends

    MPI_Request *send_reqs = (MPI_Request *) calloc(nr_send_ranks, sizeof(*send_reqs));
    int nr_send_reqs = 0;

    int **nr_particles_send_by_ri = (int **) calloc(nr_send_ranks, sizeof(*nr_particles_send_by_ri));
    for (int ri = 0; ri < nr_send_ranks; ri++) {
      struct by_ri *send = &send_by_ri[ri];
      nr_particles_send_by_ri[ri] = (int *) calloc(send->nr_patches, sizeof(*nr_particles_send_by_ri[ri]));

      for (int pi = 0; pi < send->nr_patches; pi++) {
	nr_particles_send_by_ri[ri][pi] = n_prts_by_patch_old[send->pi_to_patch[pi]];
      }

      if (send->rank != mpi_rank) {
	//mprintf("send -> %d (len %d)\n", r, nr_patches_send_by_ri[ri]);
	MPI_Isend(nr_particles_send_by_ri[ri], send->nr_patches, MPI_INT,
		  send->rank, 10, comm, &send_reqs[nr_send_reqs++]);
      }
    }
    assert(nr_send_reqs <= nr_send_ranks);

    // copy local particle numbers
    {
      int send_ri = send_rank_to_ri[mpi_rank];
      int recv_ri = recv_rank_to_ri[mpi_rank];
      if (send_ri < 0) { // no local patches to copy
	assert(recv_ri < 0);
      } else {
	assert(send_by_ri[send_ri].nr_patches == recv_by_ri[recv_ri].nr_patches);
	for (int n = 0; n < send_by_ri[send_ri].nr_patches; n++) {
	  nr_particles_recv_by_ri[recv_ri][n] = nr_particles_send_by_ri[send_ri][n];
	}
      }
    }

    MPI_Waitall(nr_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);

    // update from received data

    for (int ri = 0; ri < nr_recv_ranks; ri++) {
      struct by_ri *recv = &recv_by_ri[ri];
      for (int pi = 0; pi < recv_by_ri[ri].nr_patches; pi++) {
	n_prts_by_patch_new[recv->pi_to_patch[pi]] = nr_particles_recv_by_ri[ri][pi];
      }
    }

    // clean up recv

    for (int ri = 0; ri < nr_recv_ranks; ri++) {
      free(nr_particles_recv_by_ri[ri]);
    }
    free(nr_particles_recv_by_ri);
    free(recv_reqs);

    MPI_Waitall(nr_send_reqs, send_reqs, MPI_STATUSES_IGNORE);

    // clean up send

    for (int ri = 0; ri < nr_send_ranks; ri++) {
      free(nr_particles_send_by_ri[ri]);
    }
    free(nr_particles_send_by_ri);
    free(send_reqs);

    // return result

    prof_stop(pr);
    return n_prts_by_patch_new;
  }
};

// ======================================================================
// Balance_

template<typename Mparticles, typename MfieldsState, typename Mfields>
struct Balance_ : BalanceBase
{
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mparticles::real_t;

  Balance_(int every, double factor_fields=1., bool print_loads=false, bool write_loads=false)
    : every_(every), factor_fields_(factor_fields),
      print_loads_(print_loads), write_loads_(write_loads)
  {}

  ~Balance_()
  {
    delete[] psc_balance_comp_time_by_patch;
  }
  
  std::vector<uint> initial(const std::vector<uint>& n_prts_by_patch_old) override
  {
    auto loads = get_loads_initial(*ggrid, n_prts_by_patch_old);
    return balance(loads, nullptr, n_prts_by_patch_old);
  }

  void operator()(MparticlesBase& mp) override
  {
    static int st_time_balance;
    if (!st_time_balance) {
      st_time_balance = psc_stats_register("time balancing");
    }

    psc_stats_start(st_time_balance);
    auto loads = get_loads(mp.grid(), mp);
    balance(loads, &mp);
    psc_stats_stop(st_time_balance);
  }

private:
  std::vector<double> get_loads_initial(const Grid_t& grid, const std::vector<uint>& n_prts_by_patch)
  {
    std::vector<double> loads;
    loads.reserve(n_prts_by_patch.size());

    const int *ldims = grid.ldims;
    for (auto n_prts : n_prts_by_patch) {
      loads.push_back(n_prts + factor_fields_ * ldims[0] * ldims[1] * ldims[2]);
    }

    return loads;
  }

  std::vector<double> get_loads(const Grid_t& grid, MparticlesBase& mp)
  {
    uint n_prts_by_patch[mp.n_patches()];
    mp.get_size_all(n_prts_by_patch);

    std::vector<double> loads;
    loads.reserve(mp.n_patches());
    for (int p = 0; p < mp.n_patches(); p++) {
      double load;
      if (factor_fields_ >= 0.) {
	const int *ldims = grid.ldims;
	load = n_prts_by_patch[p] + factor_fields_ * ldims[0] * ldims[1] * ldims[2];
	//mprintf("loads p %d %g %g ratio %g\n", p, loads[p], comp_time, loads[p] / comp_time);
      } else {
	load = psc_balance_comp_time_by_patch[p];
      }
#if 0
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      load = 1 + rank;
#endif
      loads.push_back(load);
    }
    return loads;
  }

  std::vector<double> gather_loads(const Grid_t& grid, std::vector<double> loads)
  {
    const MrcDomain& domain = grid.mrc_domain_;
    MPI_Comm comm = domain.comm();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // gather nr_patches for all procs on proc 0
    int *nr_patches_all = NULL;
    if (rank == 0) {
      nr_patches_all = (int *) calloc(size, sizeof(*nr_patches_all));
    }
    int n_patches = loads.size();
    MPI_Gather(&n_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);

    // gather loads for all patches on proc 0
    int *displs = NULL;
    std::vector<double> loads_all;
    if (rank == 0) {
      displs = (int *) calloc(size, sizeof(*displs));
      int off = 0;
      for (int i = 0; i < size; i++) {
	displs[i] = off;
	off += nr_patches_all[i];
      }
      int n_global_patches = domain.nGlobalPatches();
      loads_all.resize(n_global_patches);
    }
    MPI_Gatherv(loads.data(), n_patches, MPI_DOUBLE, loads_all.data(), nr_patches_all, displs,
		MPI_DOUBLE, 0, comm);

    if (rank == 0) {
#if 0
      int rl = 0;
      char s[20]; sprintf(s, "loads-%06d.asc", grid.timestep());
      FILE *f = fopen(s, "w");
      for (int p = 0; p < *p_nr_global_patches; p++) {
	if (rl < size - 1 && p >= displs[rl+1]) {
	  rl++;
	}
	fprintf(f, "%d %g %d\n", p, loads_all[p], rl);
      }
      fclose(f);
#endif
      free(nr_patches_all);
      free(displs);
    }

    return loads_all;
  }

  int find_best_mapping(const Grid_t& grid, const std::vector<double>& loads_all)
  {
    const MrcDomain& domain = grid.mrc_domain_;
    
    MPI_Comm comm = domain.comm();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int nr_patches_old = domain.nPatches();

    std::vector<int> nr_patches_all_old(size);
    MPI_Allgather(&nr_patches_old, 1, MPI_INT, nr_patches_all_old.data(), 1, MPI_INT, comm);

    std::vector<int> nr_patches_all_new(size);

    if (rank == 0) { // do the mapping on proc 0
      std::vector<double> capability(size);
      for (int p = 0; p < size; p++) {
	capability[p] = capability_default(p);
      }
      double loads_sum = 0.;
      for (int i = 0; i < loads_all.size(); i++) {
	loads_sum += loads_all[i];
      }
      double capability_sum = 0.;
      for (int p = 0; p < size; p++) {
	capability_sum += capability[p];
      }
      double load_target = loads_sum / capability_sum;
      mprintf("psc_balance: loads_sum %g capability_sum %g load_target %g\n",
	      loads_sum, capability_sum, load_target);
      int p = 0, nr_new_patches = 0;
      double load = 0.;
      double next_target = load_target * capability[0];
      for (int i = 0; i < loads_all.size(); i++) {
	load += loads_all[i];
	nr_new_patches++;
	if (p < size - 1) {
	  // if load limit is reached, or we have only as many patches as processors left
	  if (load > next_target || size - p >= loads_all.size() - i) {
	    double above_target = load - next_target;
	    double below_target = next_target - (load - loads_all[i]);
	    if (above_target > below_target && nr_new_patches > 1) {
	      nr_patches_all_new[p] = nr_new_patches - 1;
	      nr_new_patches = 1;
	    } else {
	      nr_patches_all_new[p] = nr_new_patches;
	      nr_new_patches = 0;
	    }
	    p++;
	    next_target += load_target * capability[p];
	  }
	}
	// last proc takes what's left
	if (i == loads_all.size() - 1) {
	  nr_patches_all_new[size - 1] = nr_new_patches;
	}
      }
    
      int pp = 0;
      double min_diff = 0, max_diff = 0;
      for (int p = 0; p < size; p++) {
	double load = 0.;
	for (int i = 0; i < nr_patches_all_new[p]; i++) {
	  load += loads_all[pp++];
	  if (print_loads_) {
	    mprintf("  pp %d load %g : %g\n", pp-1, loads_all[pp-1], load);
	  }
	}
	double diff = load - load_target * capability[p];
	if (print_loads_) {
	  mprintf("p %d # = %d load %g / %g : diff %g %%\n", p, nr_patches_all_new[p],
		  load, load_target * capability[p], 100. * diff / (load_target * capability[p]));
	}
	if (diff < min_diff) {
	  min_diff = diff;
	}
	if (diff > max_diff) {
	  max_diff = diff;
	}
      }
      mprintf("psc_balance: achieved target %g (%g %% -- %g %%)\n", load_target,
	      100 * min_diff / load_target, 100 * max_diff / load_target);

      if (write_loads_) {
	int gp = 0;
	char s[20]; sprintf(s, "loads2-%06d.asc", grid.timestep());
	FILE *f = fopen(s, "w");
	for (int r = 0; r < size; r++) {
	  for (int p = 0; p < nr_patches_all_new[r]; p++) {
	    fprintf(f, "%d %g %d\n", gp, loads_all[gp], r);
	    gp++;
	  }
	}
	fclose(f);
      }

      if (nr_patches_all_new == nr_patches_all_old) {
	std::fill(nr_patches_all_new.begin(), nr_patches_all_new.end(), -1); // unchanged mapping, no communication etc needed
      }
    }

    // then scatter
    int nr_patches_new;
    MPI_Scatter(nr_patches_all_new.data(), 1, MPI_INT, &nr_patches_new, 1, MPI_INT,
		0, comm);
    return nr_patches_new;
  }

  void communicate_particles(struct communicate_ctx *ctx, Mparticles& mp_old, Mparticles& mp_new,
			     const std::vector<uint>& n_prts_by_patch_new)
  {
    static int pr, pr_A, pr_B, pr_C, pr_D;
    if (!pr) {
      pr   = prof_register("comm prts", 1., 0, 0);
      pr_A = prof_register("comm prts A", 1., 0, 0);
      pr_B = prof_register("comm prts B", 1., 0, 0);
      pr_C = prof_register("comm prts C", 1., 0, 0);
      pr_D = prof_register("comm prts D", 1., 0, 0);
    }

    prof_start(pr);

    prof_start(pr_A);
    mp_new.reserve_all(n_prts_by_patch_new.data());
    mp_new.resize_all(n_prts_by_patch_new.data());

    assert(sizeof(particle_t) % sizeof(real_t) == 0); // FIXME

    MPI_Datatype mpi_dtype = Mparticles_traits<Mparticles>::mpi_dtype();
    // recv for new local patches
    MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new]();
    int nr_recv_reqs = 0;

    for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
      struct by_ri *recv = &ctx->recv_by_ri[ri];
      if (recv->rank == ctx->mpi_rank) {
	continue;
      }

      for (int pi = 0; pi < recv->nr_patches; pi++) {
	int p = recv->pi_to_patch[pi];
	auto& prts_new = mp_new[p];
	int nn = prts_new.size() * (sizeof(particle_t)  / sizeof(real_t));
	MPI_Irecv(&*prts_new.begin(), nn, mpi_dtype, recv->rank,
		  pi, ctx->comm, &recv_reqs[nr_recv_reqs++]);
      }
    }
    prof_stop(pr_A);

    prof_start(pr_B);
    // send from old local patches
    MPI_Request *send_reqs = new MPI_Request[ctx->nr_patches_old]();
    int nr_send_reqs = 0;

    for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
      struct by_ri *send = &ctx->send_by_ri[ri];
      if (send->rank == ctx->mpi_rank) {
	continue;
      }

      for (int pi = 0; pi < send->nr_patches; pi++) {
	int p = send->pi_to_patch[pi];
	auto& prts_old = mp_old[p];
	int nn = prts_old.size() * (sizeof(particle_t)  / sizeof(real_t));
	//mprintf("A send -> %d tag %d (patch %d)\n", send->rank, pi, p);
	MPI_Isend(&*prts_old.begin(), nn, mpi_dtype, send->rank,
		  pi, ctx->comm, &send_reqs[nr_send_reqs++]);
      }
    }
    prof_stop(pr_B);

    prof_start(pr_C);
    // local particles
    // OPT: could keep the alloced arrays, just move pointers...
    for (int p = 0; p < ctx->nr_patches_new; p++) {
      if (ctx->recv_info[p].rank != ctx->mpi_rank) {
	continue;
      }

      auto& prts_old = mp_old[ctx->recv_info[p].patch];
      auto& prts_new = mp_new[p];
      assert(prts_old.size() == prts_new.size());
#if 1
      for (int n = 0; n < prts_new.size(); n++) {
	prts_new[n] = prts_old[n];
      }
#else
      // FIXME, this needs at least a proper interface -- if not separately alloc'ed, bad things
      // are going to happen
      free(c_new->particles);
      c_new->particles = c_old->particles;
      c_new->n_alloced = c_old->n_alloced;
      c_old->particles = NULL; // prevent from being freed
      c_old->n_alloced = 0;
#endif
    }
    prof_stop(pr_C);
  
    prof_start(pr_D);
    MPI_Waitall(nr_send_reqs, send_reqs, MPI_STATUSES_IGNORE);
    MPI_Waitall(nr_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);
    delete[] send_reqs;
    delete[] recv_reqs;
    prof_stop(pr_D);

    prof_stop(pr);
  }

  void communicate_fields(struct communicate_ctx *ctx, Mfields& mf_old, Mfields& mf_new)
  {
    MPI_Datatype mpi_dtype = Mfields_traits<Mfields>::mpi_dtype();

    // send from old local patches
    MPI_Request *send_reqs = new MPI_Request[ctx->nr_patches_old]();
    int *nr_patches_new_by_rank = new int[ctx->mpi_size]();
    for (int p = 0; p < ctx->nr_patches_old; p++) {
      int new_rank = ctx->send_info[p].rank;
      if (new_rank == ctx->mpi_rank || new_rank < 0) {
	send_reqs[p] = MPI_REQUEST_NULL;
      } else {
	fields_t flds_old = mf_old[p];
	Fields F_old(flds_old);
	int nn = flds_old.size();
	Int3 ib = flds_old.ib();
	void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
	int tag = nr_patches_new_by_rank[new_rank]++;
	MPI_Isend(addr_old, nn, mpi_dtype, new_rank, tag, ctx->comm, &send_reqs[p]);
      }
    }
    delete[] nr_patches_new_by_rank;

    // recv for new local patches
    MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new]();
    int *nr_patches_old_by_rank = new int[ctx->mpi_size]();
    for (int p = 0; p < ctx->nr_patches_new; p++) {
      int old_rank = ctx->recv_info[p].rank;
      if (old_rank == ctx->mpi_rank) {
	recv_reqs[p] = MPI_REQUEST_NULL;
      } else if (old_rank < 0) { //this patch did not exist before
	recv_reqs[p] = MPI_REQUEST_NULL;
	//Seed new data
      } else {
	fields_t flds_new = mf_new[p];
	Fields F_new(flds_new);
	int nn = flds_new.size();
	Int3 ib = flds_new.ib();
	void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
	int tag = nr_patches_old_by_rank[old_rank]++;
	MPI_Irecv(addr_new, nn, mpi_dtype, old_rank,
		  tag, ctx->comm, &recv_reqs[p]);
      }
    }
    delete[] nr_patches_old_by_rank;

    static int pr;
    if (!pr) {
      pr = prof_register("bal flds local", 1., 0, 0);
    }

    prof_start(pr);
    // local fields
    // OPT: could keep the alloced arrays, just move pointers...
    for (int p = 0; p < ctx->nr_patches_new; p++) {
      if (ctx->recv_info[p].rank != ctx->mpi_rank) {
	continue;
      }

      fields_t flds_old = mf_old[ctx->recv_info[p].patch];
      fields_t flds_new = mf_new[p];
      Fields F_old(flds_old), F_new(flds_new);
      assert(flds_old.n_comps() == flds_new.n_comps());
      assert(flds_old.size() == flds_new.size());
      int size = flds_old.size();
      Int3 ib = flds_new.ib();
      void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
      void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
      memcpy(addr_new, addr_old, size * sizeof(typename fields_t::real_t));
    }
    prof_stop(pr);

    MPI_Waitall(ctx->nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
    MPI_Waitall(ctx->nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
    delete[] send_reqs;
    delete[] recv_reqs;
  }

  void balance_particles(communicate_ctx& ctx, const Grid_t& new_grid, MparticlesBase& mp_base)
  {
    int n_patches = mp_base.n_patches();
    std::vector<uint> n_prts_by_patch_old(n_patches);
    mp_base.get_size_all(n_prts_by_patch_old.data());
    auto n_prts_by_patch_new = ctx.new_n_prts(n_prts_by_patch_old);

    if (typeid(mp_base) != typeid(Mparticles)) {
      auto& mp_old = *new Mparticles{mp_base.grid()};
      MparticlesBase::convert(mp_base, mp_old);
      mp_base.reset(new_grid); // frees memory here already
      
      auto mp_new = Mparticles{new_grid};
      communicate_particles(&ctx, mp_old, mp_new, n_prts_by_patch_new);
      delete &mp_old;
      
      MparticlesBase::convert(mp_new, mp_base);
    } else {
      auto mp_new = Mparticles{new_grid};
      auto& mp_old = dynamic_cast<Mparticles&>(mp_base);

      communicate_particles(&ctx, mp_old, mp_new, n_prts_by_patch_new);

      mp_old = std::move(mp_new);
    }
  }
  
  void balance_field(communicate_ctx& ctx, const Grid_t& new_grid, MfieldsBase& mf_base)
  {
    if (typeid(mf_base) != typeid(Mfields)) {
      auto& mf_old = *new Mfields{mf_base.grid(), mf_base.n_comps(), mf_base.ibn()};
      MfieldsBase::convert(mf_base, mf_old, 0, mf_old.n_comps());
      mf_base.reset(new_grid); // free old memory

      auto mf_new = Mfields{new_grid, mf_base.n_comps(), mf_base.ibn()};
      communicate_fields(&ctx, mf_old, mf_new);
      delete &mf_old; // delete as early as possible
      
      MfieldsBase::convert(mf_new, mf_base, 0, mf_base.n_comps());
    } else {
      auto mf_new = Mfields{new_grid, mf_base.n_comps(), mf_base.ibn()};
      auto &mf_old = dynamic_cast<Mfields&>(mf_base);

      communicate_fields(&ctx, mf_old, mf_new);

      mf_old = std::move(mf_new);
    }
  }
  
  void balance_state_field(communicate_ctx& ctx, const Grid_t& new_grid, MfieldsStateBase& mf_base)
  {
    if (typeid(mf_base) != typeid(MfieldsState)) {
      assert(0);
#if 0
      auto& mf_old = *new MfieldsState{mf_base.grid()};
      MfieldsStateBase::convert(mf_base, mf_old, 0, mf_old.n_comps());
      mf_base.reset(new_grid); // free old memory

      auto mf_new = MfieldsState{new_grid};
      communicate_fields(&ctx, mf_old.mflds(), mf_new.mflds());
      delete &mf_old; // delete as early as possible
      
      MfieldsBase::convert(mf_new.mflds(), mf_base.mflds(), 0, mf_base.n_comps());
#endif
    } else {
      // for now MfieldsStateFromMfields will have its underlying Mfields rebalanced, anyway, so all we need ot do
      // is reset the grid.
      // FIXME, that his however broken if using MfieldsState that isn't MfieldsStateFromMfields...
      mf_base.reset(new_grid);
#if 0
      auto mf_new = MfieldsState{new_grid};
      auto &mf_old = dynamic_cast<MfieldsState&>(mf_base);

      communicate_fields(&ctx, mf_old.mflds(), mf_new.mflds());

      mf_old = std::move(mf_new);
#endif
    }
  }
  
  std::vector<uint> balance(std::vector<double> loads, MparticlesBase *mp,
			    std::vector<uint> n_prts_by_patch_old = {})
  {
    static int pr_bal_load, pr_bal_ctx, pr_bal_prts, pr_bal_flds;
    if (!pr_bal_load) {
      pr_bal_load = prof_register("bal load", 1., 0, 0);
      pr_bal_ctx  = prof_register("bal ctx", 1., 0, 0);
      pr_bal_prts = prof_register("bal prts", 1., 0, 0);
      pr_bal_flds = prof_register("bal flds", 1., 0, 0);
    }

    prof_start(pr_bal_load);
    auto old_grid = ggrid;
    
    auto loads_all = gather_loads(*old_grid, loads);
    int n_patches_new = find_best_mapping(*old_grid, loads_all);

    auto new_grid = new Grid_t{old_grid->domain, old_grid->bc, old_grid->kinds,
			       old_grid->norm, old_grid->dt, n_patches_new};
    new_grid->ibn = old_grid->ibn; // FIXME, sucky ibn handling...
    new_grid->timestep_ = old_grid->timestep_;
    
    prof_stop(pr_bal_load);
    if (n_patches_new < 0) { // unchanged mapping, nothing tbd
      mpi_printf(old_grid->comm(), "***** Balance: decomposition unchanged\n");
      return n_prts_by_patch_old;
    }

    mpi_printf(old_grid->comm(), "***** Balance: new decomposition: balancing\n");
    delete[] psc_balance_comp_time_by_patch;
    psc_balance_comp_time_by_patch = new double[new_grid->n_patches()];
    
    prof_start(pr_bal_ctx);
    communicate_ctx ctx(old_grid->mrc_domain(), new_grid->mrc_domain());
    prof_stop(pr_bal_ctx);

    // particles
    std::vector<uint> n_prts_by_patch_new;
    if (mp) {
      prof_start(pr_bal_prts);
      balance_particles(ctx, *new_grid, *mp);
      prof_stop(pr_bal_prts);
    } else {
      n_prts_by_patch_new = ctx.new_n_prts(n_prts_by_patch_old);
    }

    prof_start(pr_bal_flds);
    // state field
    for (auto mf : MfieldsStateBase::instances) {
      balance_state_field(ctx, *new_grid, *mf);
    }
    
    // fields
    for (auto mf : MfieldsBase::instances) {
      balance_field(ctx, *new_grid, *mf);
    }
    prof_stop(pr_bal_flds);

    // update psc etc
    delete ggrid;
    ggrid = new_grid;
    psc_balance_generation_cnt++;

    return n_prts_by_patch_new;
  }

private:
  int every_;
  double factor_fields_;
  bool print_loads_;
  bool write_loads_;
};
