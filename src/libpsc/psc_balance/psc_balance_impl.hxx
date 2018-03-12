
#include "fields.hxx"
#include "balance.hxx"
#include "bnd_particles.hxx"
#include "bnd.hxx"

#include <mrc_profile.h>
#include <string.h>

extern double *psc_balance_comp_time_by_patch;
extern bool psc_output_fields_check_bnd;

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

static void
psc_get_loads_initial(struct psc *psc, double *loads, uint *nr_particles_by_patch)
{
  psc_foreach_patch(psc, p) {
    const int *ldims = psc->grid().ldims;
    loads[p] = nr_particles_by_patch[p] + 
      psc->balance->factor_fields * ldims[0] * ldims[1] * ldims[2];
  }
}

static void
psc_get_loads(struct psc *psc, double *loads)
{
  struct psc_mparticles *mprts = psc->particles;
  PscMparticlesBase mp(mprts);
  
  uint n_prts_by_patch[mp->n_patches()];
  mp->get_size_all(n_prts_by_patch);
  psc_foreach_patch(psc, p) {
    if (psc->balance->factor_fields >= 0.) {
      const int *ldims = psc->grid().ldims;
      loads[p] = n_prts_by_patch[p] +
	psc->balance->factor_fields * ldims[0] * ldims[1] * ldims[2];
      //mprintf("loads p %d %g %g ratio %g\n", p, loads[p], comp_time, loads[p] / comp_time);
    } else {
      double comp_time = psc_balance_comp_time_by_patch[p];
      loads[p] = comp_time;
#if 0
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      loads[p] = 1 + rank;
#endif
    }
  }
}

static int
find_best_mapping(struct psc_balance *bal, struct mrc_domain *domain,
		  int nr_global_patches, double *loads_all)
{
  MPI_Comm comm = mrc_domain_comm(domain);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int *nr_patches_all_new = NULL;

  if (rank == 0) { // do the mapping on proc 0
    double *capability = (double *) calloc(size, sizeof(*capability));
    for (int p = 0; p < size; p++) {
      capability[p] = capability_default(p);
    }
    nr_patches_all_new = (int *) calloc(size, sizeof(*nr_patches_all_new));
    double loads_sum = 0.;
    for (int i = 0; i < nr_global_patches; i++) {
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
    for (int i = 0; i < nr_global_patches; i++) {
      load += loads_all[i];
      nr_new_patches++;
      if (p < size - 1) {
	// if load limit is reached, or we have only as many patches as processors left
	if (load > next_target || size - p >= nr_global_patches - i) {
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
      if (i == nr_global_patches - 1) {
	nr_patches_all_new[size - 1] = nr_new_patches;
      }
    }
    
    int pp = 0;
    double min_diff = 0, max_diff = 0;
    for (int p = 0; p < size; p++) {
      double load = 0.;
      for (int i = 0; i < nr_patches_all_new[p]; i++) {
	load += loads_all[pp++];
	if (bal->print_loads) {
	  mprintf("  pp %d load %g : %g\n", pp-1, loads_all[pp-1], load);
	}
      }
      double diff = load - load_target * capability[p];
      if (bal->print_loads) {
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

    if (bal->write_loads) {
      int gp = 0;
      char s[20]; sprintf(s, "loads2-%06d.asc", ppsc->timestep);
      FILE *f = fopen(s, "w");
      for (int r = 0; r < size; r++) {
	for (int p = 0; p < nr_patches_all_new[r]; p++) {
	  fprintf(f, "%d %g %d\n", gp, loads_all[gp], r);
	  gp++;
	}
      }
      fclose(f);
    }

    free(capability);
  }
  // then scatter
  int nr_patches_new;
  MPI_Scatter(nr_patches_all_new, 1, MPI_INT, &nr_patches_new, 1, MPI_INT,
	      0, comm);
  free(nr_patches_all_new);
  return nr_patches_new;
}

static double *
gather_loads(struct mrc_domain *domain, double *loads, int nr_patches,
	     int *p_nr_global_patches)
{
  MPI_Comm comm = mrc_domain_comm(domain);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // gather nr_patches for all procs on proc 0
  int *nr_patches_all = NULL;
  if (rank == 0) {
    nr_patches_all = (int *) calloc(size, sizeof(*nr_patches_all));
  }
  MPI_Gather(&nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);

  // gather loads for all patches on proc 0
  int *displs = NULL;
  double *loads_all = NULL;
  if (rank == 0) {
    displs = (int *) calloc(size, sizeof(*displs));
    int off = 0;
    for (int i = 0; i < size; i++) {
      displs[i] = off;
      off += nr_patches_all[i];
    }
    mrc_domain_get_nr_global_patches(domain, p_nr_global_patches);
	  
    loads_all = (double *) calloc(*p_nr_global_patches, sizeof(*loads_all));
  }
  MPI_Gatherv(loads, nr_patches, MPI_DOUBLE, loads_all, nr_patches_all, displs,
	      MPI_DOUBLE, 0, comm);

  if (rank == 0) {
#if 0
    int rl = 0;
    char s[20]; sprintf(s, "loads-%06d.asc", ppsc->timestep);
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

static void
communicate_setup(struct communicate_ctx *ctx, struct mrc_domain *domain_old,
		  struct mrc_domain *domain_new)
{
  ctx->comm = mrc_domain_comm(domain_old);
  MPI_Comm_rank(ctx->comm, &ctx->mpi_rank);
  MPI_Comm_size(ctx->comm, &ctx->mpi_size);

  mrc_domain_get_patches(domain_old, &ctx->nr_patches_old);
  mrc_domain_get_patches(domain_new, &ctx->nr_patches_new);

  ctx->send_info = (struct send_info *) calloc(ctx->nr_patches_old, sizeof(*ctx->send_info));
  ctx->recv_info = (struct recv_info *) calloc(ctx->nr_patches_new, sizeof(*ctx->recv_info));

  for (int p = 0; p < ctx->nr_patches_old; p++) {
    struct mrc_patch_info info_old, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info_old);
    mrc_domain_get_level_idx3_patch_info(domain_new, info_old.level, info_old.idx3, &info_new);
    ctx->send_info[p].rank  = info_new.rank;
    ctx->send_info[p].patch = info_new.patch;
  }
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    struct mrc_patch_info info_old, info_new;
    mrc_domain_get_local_patch_info(domain_new, p, &info_new);
    mrc_domain_get_level_idx3_patch_info(domain_old, info_new.level, info_new.idx3, &info_old);
    ctx->recv_info[p].rank  = info_old.rank;
    ctx->recv_info[p].patch = info_old.patch;
  }

  // maps rank <-> rank index

  ctx->send_rank_to_ri = (int *) malloc(ctx->mpi_size * sizeof(*ctx->send_rank_to_ri));
  ctx->recv_rank_to_ri = (int *) malloc(ctx->mpi_size * sizeof(*ctx->recv_rank_to_ri));
  for (int r = 0; r < ctx->mpi_size; r++) {
    ctx->send_rank_to_ri[r] = -1;
    ctx->recv_rank_to_ri[r] = -1;
  }

  ctx->nr_send_ranks = 0;
  ctx->nr_recv_ranks = 0;
  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int send_rank = ctx->send_info[p].rank;
    if (send_rank >= 0) {
      if (ctx->send_rank_to_ri[send_rank] < 0) {
	ctx->send_rank_to_ri[send_rank] = ctx->nr_send_ranks++;
      }
    }
  }
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int recv_rank = ctx->recv_info[p].rank;
    if (recv_rank >= 0) {
      if (ctx->recv_rank_to_ri[recv_rank] < 0) {
	ctx->recv_rank_to_ri[recv_rank] = ctx->nr_recv_ranks++;
      }
    }
  }

  ctx->send_by_ri = (struct by_ri *) calloc(ctx->nr_send_ranks, sizeof(*ctx->send_by_ri));
  ctx->recv_by_ri = (struct by_ri *) calloc(ctx->nr_recv_ranks, sizeof(*ctx->recv_by_ri));

  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int send_rank = ctx->send_info[p].rank;
    if (send_rank >= 0) {
      ctx->send_by_ri[ctx->send_rank_to_ri[send_rank]].rank = send_rank;
    }
  }
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int recv_rank = ctx->recv_info[p].rank;
    if (recv_rank >= 0) {
      ctx->recv_by_ri[ctx->recv_rank_to_ri[recv_rank]].rank = recv_rank;
    }
  }

  /* for (int ri = 0; ri < ctx->nr_send_ranks; ri++) { */
  /*   mprintf("send -> %d (%d)\n", ctx->send_by_ri[ri].rank, ri); */
  /* } */

  // count number of patches sent to each rank

  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int send_rank = ctx->send_info[p].rank;
    if (send_rank >= 0) {
      ctx->send_by_ri[ctx->send_rank_to_ri[send_rank]].nr_patches++;
    }
  }

  // map send patch index by ri back to local patch number

  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    ctx->send_by_ri[ri].pi_to_patch = (int *) calloc(ctx->send_by_ri[ri].nr_patches, sizeof(*ctx->send_by_ri[ri].pi_to_patch));
    ctx->send_by_ri[ri].nr_patches = 0;
  }

  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int send_rank = ctx->send_info[p].rank;
    if (send_rank < 0) {
      continue;
    }
    int ri = ctx->send_rank_to_ri[send_rank];
    int pi = ctx->send_by_ri[ri].nr_patches++;
    ctx->send_by_ri[ri].pi_to_patch[pi] = p;
  }

  // count number of patches received from each rank

  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int recv_rank = ctx->recv_info[p].rank;
    if (recv_rank >= 0) {
      int ri = ctx->recv_rank_to_ri[recv_rank];
      ctx->recv_by_ri[ri].nr_patches++;
    }
  }

  // map received patch index by ri back to local patch number

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    ctx->recv_by_ri[ri].pi_to_patch = (int *) calloc(ctx->recv_by_ri[ri].nr_patches, sizeof(*ctx->recv_by_ri[ri].pi_to_patch));
    ctx->recv_by_ri[ri].nr_patches = 0;
  }

  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int recv_rank = ctx->recv_info[p].rank;
    if (recv_rank < 0) {
      continue;
    }
    int ri = ctx->recv_rank_to_ri[recv_rank];
    int pi = ctx->recv_by_ri[ri].nr_patches++;
    ctx->recv_by_ri[ri].pi_to_patch[pi] = p;
  }

}

static void
communicate_free(struct communicate_ctx *ctx)
{
  free(ctx->send_info);
  free(ctx->recv_info);

  free(ctx->send_rank_to_ri);
  free(ctx->recv_rank_to_ri);

  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    free(ctx->send_by_ri[ri].pi_to_patch);
  }
  free(ctx->send_by_ri);

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    free(ctx->recv_by_ri[ri].pi_to_patch);
  }
  free(ctx->recv_by_ri);
}

static void
communicate_new_nr_particles(struct communicate_ctx *ctx, uint **p_nr_particles_by_patch)
{
  static int pr;
  if (!pr) {
    pr   = prof_register("comm nr prts", 1., 0, 0);
  }

  prof_start(pr);

  uint *nr_particles_by_patch_old = *p_nr_particles_by_patch;
  uint *nr_particles_by_patch_new = (uint *) calloc(ctx->nr_patches_new,
						   sizeof(nr_particles_by_patch_new));
  // post receives 

  MPI_Request *recv_reqs = (MPI_Request *) calloc(ctx->nr_patches_new, sizeof(*recv_reqs));
  int nr_recv_reqs = 0;

  int **nr_particles_recv_by_ri = (int **) calloc(ctx->nr_recv_ranks, sizeof(*nr_particles_recv_by_ri));
  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    struct by_ri *recv = &ctx->recv_by_ri[ri];
    nr_particles_recv_by_ri[ri] = (int *) calloc(recv->nr_patches, sizeof(*nr_particles_recv_by_ri[ri]));
    
    if (recv->rank != ctx->mpi_rank) {
      //mprintf("recv <- %d (len %d)\n", r, nr_patches_recv_by_ri[ri]);
      MPI_Irecv(nr_particles_recv_by_ri[ri], recv->nr_patches, MPI_INT,
		recv->rank, 10, ctx->comm, &recv_reqs[nr_recv_reqs++]);
    }
  }

  // post sends

  MPI_Request *send_reqs = (MPI_Request *) calloc(ctx->nr_send_ranks, sizeof(*send_reqs));
  int nr_send_reqs = 0;

  int **nr_particles_send_by_ri = (int **) calloc(ctx->nr_send_ranks, sizeof(*nr_particles_send_by_ri));
  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    struct by_ri *send = &ctx->send_by_ri[ri];
    nr_particles_send_by_ri[ri] = (int *) calloc(send->nr_patches, sizeof(*nr_particles_send_by_ri[ri]));

    for (int pi = 0; pi < send->nr_patches; pi++) {
      nr_particles_send_by_ri[ri][pi] = nr_particles_by_patch_old[send->pi_to_patch[pi]];
    }

    if (send->rank != ctx->mpi_rank) {
      //mprintf("send -> %d (len %d)\n", r, nr_patches_send_by_ri[ri]);
      MPI_Isend(nr_particles_send_by_ri[ri], send->nr_patches, MPI_INT,
		send->rank, 10, ctx->comm, &send_reqs[nr_send_reqs++]);
    }
  }
  assert(nr_send_reqs <= ctx->nr_send_ranks);

  // copy local particle numbers
  {
    int send_ri = ctx->send_rank_to_ri[ctx->mpi_rank];
    int recv_ri = ctx->recv_rank_to_ri[ctx->mpi_rank];
    if (send_ri < 0) { // no local patches to copy
      assert(recv_ri < 0);
    } else {
      assert(ctx->send_by_ri[send_ri].nr_patches == ctx->recv_by_ri[recv_ri].nr_patches);
      for (int n = 0; n < ctx->send_by_ri[send_ri].nr_patches; n++) {
	nr_particles_recv_by_ri[recv_ri][n] = nr_particles_send_by_ri[send_ri][n];
      }
    }
  }

  MPI_Waitall(nr_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);

  // update from received data

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    struct by_ri *recv = &ctx->recv_by_ri[ri];
    for (int pi = 0; pi < ctx->recv_by_ri[ri].nr_patches; pi++) {
      nr_particles_by_patch_new[recv->pi_to_patch[pi]] = nr_particles_recv_by_ri[ri][pi];
    }
  }

  // clean up recv

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    free(nr_particles_recv_by_ri[ri]);
  }
  free(nr_particles_recv_by_ri);
  free(recv_reqs);

  MPI_Waitall(nr_send_reqs, send_reqs, MPI_STATUSES_IGNORE);

  // clean up send

  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    free(nr_particles_send_by_ri[ri]);
  }
  free(nr_particles_send_by_ri);
  free(send_reqs);

  // return result

  free(*p_nr_particles_by_patch);
  *p_nr_particles_by_patch = nr_particles_by_patch_new;

  prof_stop(pr);
}

template<typename MP, typename MF>
struct Balance_ : BalanceBase
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;
  using particle_t = typename mparticles_t::particle_t;
  using real_t = typename mparticles_t::real_t;
  using Mparticles = typename mparticles_t::sub_t;
  using Mfields = typename mfields_t::sub_t;

  void communicate_particles(struct psc_balance *bal, struct communicate_ctx *ctx,
			     mparticles_t mprts_old, mparticles_t mprts_new,
			     uint *nr_particles_by_patch_new)
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
    // FIXME, use _all
    for (int p = 0; p < mprts_new->n_patches(); p++) {
      mprts_new[p].reserve(nr_particles_by_patch_new[p]);
      mprts_new[p].resize(nr_particles_by_patch_new[p]);
    }

    assert(sizeof(particle_t) % sizeof(real_t) == 0); // FIXME

    MPI_Datatype mpi_dtype = mparticles_traits<mparticles_t>::mpi_dtype();
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
	auto& prts_new = mprts_new[p];
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
	auto& prts_old = mprts_old[p];
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

      auto& prts_old = mprts_old[ctx->recv_info[p].patch];
      auto& prts_new = mprts_new[p];
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

  void communicate_fields(struct psc_balance *bal, struct communicate_ctx *ctx,
			  mfields_t mf_old, mfields_t mf_new)
  {
    //HACK: Don't communicate output fields if they don't correspond to the domain
    //This is needed e.g. for the boosted output which handles its MPI communication internally
    //printf("Field: %s\n", flds->f[0].name);

    if (ctx->nr_patches_old != mf_old->n_patches()) return;
	
    assert(ctx->nr_patches_old == mf_old->n_patches());
    assert(ctx->nr_patches_old > 0);
  
    MPI_Datatype mpi_dtype = fields_traits<fields_t>::mpi_dtype();

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

  void balance_field(struct psc_balance *bal, struct communicate_ctx* ctx,
		     struct psc* psc, struct mrc_domain *domain_new, const Grid_t *new_grid,
		     struct psc_mfields **p_mflds)
  {
    PscMfieldsBase mflds_base_old{*p_mflds};

    struct psc_mfields *flds_base_new;
    flds_base_new = psc_mfields_create(mrc_domain_comm(domain_new));
    psc_mfields_set_type(flds_base_new, psc_mfields_type(mflds_base_old.mflds()));
    psc_mfields_set_name(flds_base_new, psc_mfields_name(mflds_base_old.mflds()));
    psc_mfields_set_param_int(flds_base_new, "nr_fields", mflds_base_old->n_comps());
    psc_mfields_set_param_int3(flds_base_new, "ibn", mflds_base_old.mflds()->ibn);
    flds_base_new->grid = new_grid;
    psc_mfields_setup(flds_base_new);
    for (int m = 0; m < mflds_base_old->n_comps(); m++) {
      const char *s = psc_mfields_comp_name(mflds_base_old.mflds(), m);
      if (s) {
	psc_mfields_set_comp_name(flds_base_new, m, s);
      }
    }
    PscMfieldsBase mflds_base_new{flds_base_new};
    
    // FIXME, need to move up to avoid keeping two copies of CUDA fields on GPU
    mfields_t mflds_old = mflds_base_old.get_as<mfields_t>(0, mflds_base_old->n_comps());
    if (mflds_old.mflds() != mflds_base_old.mflds()) { 
      mflds_base_old->~MfieldsBase();
    }
    
    mfields_t mflds_new = mflds_base_new.get_as<mfields_t>(0, 0);
    communicate_fields(bal, ctx, mflds_old, mflds_new);
    mflds_new.put_as(mflds_base_new, 0, mflds_base_new->n_comps());
    
    if (mflds_old.mflds() == mflds_base_old.mflds()) {
      mflds_old.put_as(mflds_base_old, 0, 0);
      // psc_mfields_destroy(mflds_base_old.mflds());
      // *p_mflds = mflds_base_new.mflds();
      mflds_base_old->~MfieldsBase();
    }
    memcpy((char*) mflds_base_old.sub(), (char*) mflds_base_new.sub(),
	   mflds_base_old.mflds()->obj.ops->size);
  }
  
  void initial(struct psc_balance *bal, struct psc *psc, uint*& n_prts_by_patch) override
  {
    struct mrc_domain *domain_old = psc->mrc_domain;

    int nr_patches;
    mrc_domain_get_patches(domain_old, &nr_patches);
    double *loads = (double *) calloc(nr_patches, sizeof(*loads));
    psc_get_loads_initial(psc, loads, n_prts_by_patch);

    int nr_global_patches;
    double *loads_all = gather_loads(domain_old, loads, nr_patches,
				     &nr_global_patches);
    free(loads);

    int nr_patches_new = find_best_mapping(bal, domain_old, nr_global_patches,
					   loads_all);

    free(loads_all);

    struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_new);
    //  mrc_domain_view(domain_new);
    delete psc->grid_;
    psc->grid_ = psc->make_grid(domain_new);
    psc_balance_comp_time_by_patch = (double *) calloc(nr_patches_new,
						       sizeof(*psc_balance_comp_time_by_patch));

    struct communicate_ctx _ctx, *ctx = &_ctx;
    communicate_setup(ctx, domain_old, domain_new);

    communicate_new_nr_particles(ctx, &n_prts_by_patch);

    // ----------------------------------------------------------------------
    // fields

    struct psc_mfields_list_entry *p;
    __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
      balance_field(bal, ctx, psc, domain_new, &psc->grid(), p->flds_p);
    }

    communicate_free(ctx);

    psc->mrc_domain = domain_new;
    auto bndp = PscBndParticlesBase(psc->bnd_particles);
    bndp.reset();
    psc_output_fields_check_bnd = true;
    mrc_domain_destroy(domain_old);
  }

  void operator()(struct psc_balance *bal, struct psc *psc) override
  {
    // FIXME, way too much duplication from the above
    if (bal->every <= 0)
      return;

    if (psc->timestep == 0 || psc->timestep % bal->every != 0)
      return;

    static int st_time_balance;
    if (!st_time_balance) {
      st_time_balance = psc_stats_register("time balancing");
    }
    static int pr_bal_gather, pr_bal_decomp_A, pr_bal_decomp_B, pr_bal_decomp_C, pr_bal_decomp_D,
      pr_bal_prts, pr_bal_flds, pr_bal_flds_comm, pr_bal_ctx, pr_bal_flds_A, pr_bal_flds_B, pr_bal_flds_C;
    static int pr_bal_prts_A, pr_bal_prts_B, pr_bal_prts_C, pr_bal_prts_B1;
    if (!pr_bal_gather) {
      pr_bal_gather = prof_register("bal gather", 1., 0, 0);
      pr_bal_decomp_A = prof_register("bal decomp A", 1., 0, 0);
      pr_bal_decomp_B = prof_register("bal decomp B", 1., 0, 0);
      pr_bal_decomp_C = prof_register("bal decomp C", 1., 0, 0);
      pr_bal_decomp_D = prof_register("bal decomp D", 1., 0, 0);
      pr_bal_ctx = prof_register("bal ctx", 1., 0, 0);
      pr_bal_prts = prof_register("bal prts", 1., 0, 0);
      pr_bal_prts_A = prof_register("bal prts A", 1., 0, 0);
      pr_bal_prts_B = prof_register("bal prts B", 1., 0, 0);
      pr_bal_prts_B1 = prof_register("bal prts B1", 1., 0, 0);
      pr_bal_prts_C = prof_register("bal prts C", 1., 0, 0);
      pr_bal_flds = prof_register("bal flds", 1., 0, 0);
      pr_bal_flds_comm = prof_register("bal flds comm", 1., 0, 0);
      pr_bal_flds_A = prof_register("bal flds A", 1., 0, 0);
      pr_bal_flds_B = prof_register("bal flds B", 1., 0, 0);
      pr_bal_flds_C = prof_register("bal flds C", 1., 0, 0);
    }

    psc_stats_start(st_time_balance);
    prof_start(pr_bal_gather);
    struct mrc_domain *domain_old = psc->mrc_domain;

    int nr_patches;
    mrc_domain_get_patches(domain_old, &nr_patches);
    double *loads = (double *) calloc(nr_patches, sizeof(*loads));
    psc_get_loads(psc, loads);

    int nr_global_patches;
    double *loads_all = gather_loads(domain_old, loads, nr_patches,
				     &nr_global_patches);
    free(loads);
    prof_stop(pr_bal_gather);

    prof_start(pr_bal_decomp_A);
    int nr_patches_new = find_best_mapping(bal, domain_old, nr_global_patches,
					   loads_all);
    prof_stop(pr_bal_decomp_A);

    prof_start(pr_bal_decomp_B);
    free(loads_all);

    struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_new);
    prof_stop(pr_bal_decomp_B);
    prof_start(pr_bal_decomp_C);
    //  mrc_domain_view(domain_new);
    Grid_t* new_grid = psc->make_grid(domain_new);
    prof_stop(pr_bal_decomp_C);
    prof_start(pr_bal_decomp_D);
    free(psc_balance_comp_time_by_patch);
    psc_balance_comp_time_by_patch = (double *) calloc(nr_patches_new,// leaked the very last time
						       sizeof(*psc_balance_comp_time_by_patch));
  
    //If there are no active patches, exit here
    int n_global_patches;
    mrc_domain_get_nr_global_patches(domain_new, &n_global_patches);
    if(n_global_patches < 1) abort();
    prof_stop(pr_bal_decomp_D);

    // OPT: if local patches didn't change at all, no need to do anything...

    prof_start(pr_bal_ctx);
    struct communicate_ctx _ctx, *ctx = &_ctx;
    communicate_setup(ctx, domain_old, domain_new);
    prof_stop(pr_bal_ctx);

    // ----------------------------------------------------------------------
    // particles

    prof_start(pr_bal_prts);

    prof_start(pr_bal_prts_A);
    PscMparticlesBase mprts_old_base(psc->particles);
    uint *nr_particles_by_patch = (uint *) calloc(nr_patches, sizeof(*nr_particles_by_patch));
    mprts_old_base->get_size_all(nr_particles_by_patch);
    prof_stop(pr_bal_prts_A);

    communicate_new_nr_particles(ctx, &nr_particles_by_patch);

    prof_start(pr_bal_prts_B);
    // alloc new particles
    mrc_obj_ops *ops = mrc_obj_find_subclass_ops(reinterpret_cast<mrc_class*>(&mrc_class_psc_mparticles),
						 psc_mparticles_type(mprts_old_base.mprts()));
    assert(ops);
						 
    struct psc_mparticles *_mprts_base_new = psc_mparticles_create(mrc_domain_comm(domain_new));
    psc_mparticles_set_type(_mprts_base_new, psc->prm.particles_base);
    _mprts_base_new->grid = new_grid;
    psc_mparticles_setup(_mprts_base_new);
    PscMparticlesBase mprts_new_base{_mprts_base_new};

    mparticles_t mprts_old = mprts_old_base.get_as<mparticles_t>();
    if (mprts_old.mprts() != mprts_old_base.mprts()) { // FIXME hacky: destroy old particles early if we just got a copy
      mprts_old_base.sub()->~MparticlesBase();
      //      psc_mparticles_destroy(psc->particles);
    }

    prof_start(pr_bal_prts_B1);
    mprts_new_base->reserve_all(nr_particles_by_patch);
    prof_stop(pr_bal_prts_B1);

    mparticles_t mprts_new = mprts_new_base.get_as<mparticles_t>(MP_DONT_COPY);
    prof_stop(pr_bal_prts_B);
    
    // communicate particles
    communicate_particles(bal, ctx, mprts_old, mprts_new, nr_particles_by_patch);

    prof_start(pr_bal_prts_C);
    free(nr_particles_by_patch);

    mprts_new.put_as(mprts_new_base, 0);

    if (mprts_old.mprts() != mprts_old_base.mprts()) {
      // can't do this because psc->particles is gone
      // psc_mparticles_put_as(mprts_old, psc->particles, MP_DONT_COPY);
      psc_mparticles_destroy(mprts_old.mprts());
    } else {
      // replace particles by redistributed ones
      // FIXME, very hacky: brutally overwrites the sub-object, maybe this could be done properly
      // with move semantics
      // psc_mparticles_destroy(psc->particles);
      // psc->particles = mprts_base_new;
      mprts_old_base.sub()->~MparticlesBase();
    }
    mprts_old_base.mprts()->grid = new_grid;
    memcpy((char*) mprts_old_base.sub(), (char*) mprts_new_base.sub(), mprts_old_base.mprts()->obj.ops->size);
    prof_stop(pr_bal_prts_C);

    prof_stop(pr_bal_prts);

    // ----------------------------------------------------------------------
    // fields

    prof_start(pr_bal_flds);
    struct psc_mfields_list_entry *p;
    __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
      balance_field(bal, ctx, psc, domain_new, new_grid, p->flds_p);
    }
    prof_stop(pr_bal_flds);
  
    communicate_free(ctx);

    delete psc->grid_;
    psc->grid_ = new_grid;
    psc->mrc_domain = domain_new;
    auto bndp = PscBndParticlesBase(psc->bnd_particles);
    bndp.reset();
    auto bnd = PscBndBase(psc->bnd);
    bnd.reset();
    psc_output_fields_check_bnd = true;
    psc_balance_generation_cnt++;
  
    mrc_domain_destroy(domain_old);

    psc_stats_stop(st_time_balance);
  }
};
