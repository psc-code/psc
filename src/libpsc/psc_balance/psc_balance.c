
#include "psc_balance_private.h"
#include "psc_bnd_fields.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

LIST_HEAD(psc_mfields_base_list);

double *psc_balance_comp_time_by_patch;

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
psc_get_loads_initial(struct psc *psc, double *loads, int *nr_particles_by_patch)
{
  psc_foreach_patch(psc, p) {
    int *ldims = psc->patch[p].ldims;
    loads[p] = nr_particles_by_patch[p] + 
      psc->balance->factor_fields * ldims[0] * ldims[1] * ldims[2];
  }
}

static void
psc_get_loads(struct psc *psc, double *loads)
{
  psc_foreach_patch(psc, p) {
    if (psc->balance->factor_fields >= 0.) {
      struct psc_particles *prts = psc_mparticles_get_patch(psc->particles, p);
      int *ldims = psc->patch[p].ldims;
      loads[p] = prts->n_part +
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
    double *capability = calloc(size, sizeof(*capability));
    for (int p = 0; p < size; p++) {
      capability[p] = capability_default(p);
    }
    nr_patches_all_new = calloc(size, sizeof(*nr_patches_all_new));
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
	mprintf("p %d # = %d load %g / %g : diff %g\n", p, nr_patches_all_new[p],
		load, load_target * capability[p], diff);
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

#define MAX(x,y) (x > y ? x : y)

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
    nr_patches_all = calloc(size, sizeof(*nr_patches_all));
  }
  MPI_Gather(&nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);

  // gather loads for all patches on proc 0
  int *displs = NULL;
  double *loads_all = NULL;
  if (rank == 0) {
    displs = calloc(size, sizeof(*displs));
    int off = 0;
    for (int i = 0; i < size; i++) {
      displs[i] = off;
      off += nr_patches_all[i];
    }
    mrc_domain_get_nr_global_patches(domain, p_nr_global_patches);
	  
    //HACK If we have a dynamic domain, assume all newly created patches have a fixed load
    if(ppsc->use_dynamic_patches) {
      int n_old_patches = *p_nr_global_patches;
      *p_nr_global_patches = bitfield3d_count_bits_set(ppsc->patchmanager.activepatches);
      
      loads_all = calloc(MAX(n_old_patches, *p_nr_global_patches), sizeof(*loads_all));
      
      for(int i = n_old_patches; i < *p_nr_global_patches; i++) {
	loads_all[i] = 1.;	//TODO Better assumption? Like take the median or sth alike...
      }
    } else {
      loads_all = calloc(*p_nr_global_patches, sizeof(*loads_all));
    }
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

void
psc_balance_communicate_particles(struct psc_balance *bal, struct communicate_ctx *ctx,
				  struct psc_mparticles *mprts_old, struct psc_mparticles *mprts_new,
				  int *nr_particles_by_patch_new)
{
  struct psc_balance_ops *ops = psc_balance_ops(bal);

  assert(ops && ops->communicate_particles);
  ops->communicate_particles(bal, ctx, mprts_old, mprts_new, nr_particles_by_patch_new);
}

static void
psc_balance_communicate_fields(struct psc_balance *bal, struct communicate_ctx *ctx,
			       struct psc_mfields *mflds_old, struct psc_mfields *mflds_new)
{
  struct psc_balance_ops *ops = psc_balance_ops(bal);

  assert(ops && ops->communicate_fields);
  ops->communicate_fields(bal, ctx, mflds_old, mflds_new);
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

  ctx->send_info = calloc(ctx->nr_patches_old, sizeof(*ctx->send_info));
  ctx->recv_info = calloc(ctx->nr_patches_new, sizeof(*ctx->recv_info));

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

  ctx->send_rank_to_ri = malloc(ctx->mpi_size * sizeof(*ctx->send_rank_to_ri));
  ctx->recv_rank_to_ri = malloc(ctx->mpi_size * sizeof(*ctx->recv_rank_to_ri));
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

  ctx->send_by_ri = calloc(ctx->nr_send_ranks, sizeof(*ctx->send_by_ri));
  ctx->recv_by_ri = calloc(ctx->nr_recv_ranks, sizeof(*ctx->recv_by_ri));

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
    ctx->send_by_ri[ri].pi_to_patch = calloc(ctx->send_by_ri[ri].nr_patches, sizeof(*ctx->send_by_ri[ri].pi_to_patch));
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
    ctx->recv_by_ri[ri].pi_to_patch = calloc(ctx->recv_by_ri[ri].nr_patches, sizeof(*ctx->recv_by_ri[ri].pi_to_patch));
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
communicate_new_nr_particles(struct communicate_ctx *ctx, int **p_nr_particles_by_patch)
{
  static int pr;
  if (!pr) {
    pr   = prof_register("comm nr prts", 1., 0, 0);
  }

  prof_start(pr);

  int *nr_particles_by_patch_old = *p_nr_particles_by_patch;
  int *nr_particles_by_patch_new = calloc(ctx->nr_patches_new,
					  sizeof(nr_particles_by_patch_new));
  // post receives 

  MPI_Request *recv_reqs = calloc(ctx->nr_patches_new, sizeof(*recv_reqs));
  int nr_recv_reqs = 0;

  int **nr_particles_recv_by_ri = calloc(ctx->nr_recv_ranks, sizeof(*nr_particles_recv_by_ri));
  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    struct by_ri *recv = &ctx->recv_by_ri[ri];
    nr_particles_recv_by_ri[ri] = calloc(recv->nr_patches, sizeof(*nr_particles_recv_by_ri[ri]));
    
    if (recv->rank != ctx->mpi_rank) {
      //mprintf("recv <- %d (len %d)\n", r, nr_patches_recv_by_ri[ri]);
      MPI_Irecv(nr_particles_recv_by_ri[ri], recv->nr_patches, MPI_INT,
		recv->rank, 10, ctx->comm, &recv_reqs[nr_recv_reqs++]);
    }
  }

  // post sends

  MPI_Request *send_reqs = calloc(ctx->nr_send_ranks, sizeof(*send_reqs));
  int nr_send_reqs = 0;

  int **nr_particles_send_by_ri = calloc(ctx->nr_send_ranks, sizeof(*nr_particles_send_by_ri));
  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    struct by_ri *send = &ctx->send_by_ri[ri];
    nr_particles_send_by_ri[ri] = calloc(send->nr_patches, sizeof(*nr_particles_send_by_ri[ri]));

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

static void psc_balance_seed_patches(struct mrc_domain *domain_old, struct mrc_domain *domain_new)
{
  int nr_patches_new;
  mrc_domain_get_patches(domain_new, &nr_patches_new);
    
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_old, info.level, info.idx3, &info_old);
    if (info_old.rank < 0)	//Patch has to be seeded
    {
      //Seed field
      double t = ppsc->timestep * ppsc->dt;
      psc_bnd_fields_setup_patch(psc_push_fields_get_bnd_fields(ppsc->push_fields), p, ppsc->flds, t);

      //Seed particles
      assert(0);
      //psc_case_init_particles_patch(ppsc->patchmanager.currentcase, p, 0);
    }
  }
}

void psc_bnd_particles_check_domain(struct psc_bnd_particles *bnd);
void psc_bnd_check_domain(struct psc_bnd *bnd);
extern bool psc_output_fields_check_bnd;

void
psc_balance_initial(struct psc_balance *bal, struct psc *psc,
		    int **p_nr_particles_by_patch)
{
  struct mrc_domain *domain_old = psc->mrc_domain;

  int nr_patches;
  mrc_domain_get_patches(domain_old, &nr_patches);
  double *loads = calloc(nr_patches, sizeof(*loads));
  psc_get_loads_initial(psc, loads, *p_nr_particles_by_patch);

  int nr_global_patches;
  double *loads_all = gather_loads(domain_old, loads, nr_patches,
				   &nr_global_patches);
  free(loads);

  int nr_patches_new = find_best_mapping(bal, domain_old, nr_global_patches,
					 loads_all);

  free(loads_all);

  free(psc->patch);
  struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_new);
  //  mrc_domain_view(domain_new);
  psc_setup_patches(psc, domain_new);
  psc_balance_comp_time_by_patch = calloc(nr_patches_new,
					  sizeof(*psc_balance_comp_time_by_patch));

  struct communicate_ctx _ctx, *ctx = &_ctx;
  communicate_setup(ctx, domain_old, domain_new);

  communicate_new_nr_particles(ctx, p_nr_particles_by_patch);

  // ----------------------------------------------------------------------
  // fields

  struct psc_balance_ops *ops = psc_balance_ops(bal);

  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
    mfields_base_t *flds_base_old = *p->flds_p;
    
    // if (flds_base_old != psc->flds) {
    //   fprintf(stderr, "WARNING: not rebalancing some extra field -- expect crash!\n");
    //   continue; // FIXME!!!
    // }
    mfields_base_t *flds_base_new;
    flds_base_new = psc_mfields_create(mrc_domain_comm(domain_new));
    psc_mfields_set_type(flds_base_new, psc_mfields_type(flds_base_old));
    psc_mfields_set_name(flds_base_new, psc_mfields_name(flds_base_old));
    psc_mfields_set_domain(flds_base_new, domain_new);
    psc_mfields_set_param_int(flds_base_new, "nr_fields", flds_base_old->nr_fields);
    psc_mfields_set_param_int(flds_base_new, "first_comp", flds_base_old->first_comp);
    psc_mfields_set_param_int3(flds_base_new, "ibn", flds_base_old->ibn);
    psc_mfields_setup(flds_base_new);
    for (int m = flds_base_old->first_comp;
	 m < flds_base_old->first_comp + flds_base_old->nr_fields; m++) {
      const char *s = psc_mfields_comp_name(flds_base_old, m);
      if (s) {
	psc_mfields_set_comp_name(flds_base_new, m, s);
      }
    }

    // FIXME, need to move up to avoid keeping two copies of CUDA fields on GPU
    struct psc_mfields *flds_old =
      psc_mfields_get_as(flds_base_old, ops->mflds_type, flds_base_old->first_comp,
			 flds_base_old->first_comp + flds_base_old->nr_fields);
    if (flds_old != flds_base_old) { 
      psc_mfields_destroy(flds_base_old);
    }

    struct psc_mfields *flds_new =
      psc_mfields_get_as(flds_base_new, ops->mflds_type, 0, 0);
    psc_balance_communicate_fields(bal, ctx, flds_old, flds_new);
    psc_mfields_put_as(flds_new, flds_base_new, flds_base_new->first_comp,
		       flds_base_new->first_comp + flds_base_new->nr_fields);

    if (flds_old == flds_base_old) {
      psc_mfields_put_as(flds_old, flds_base_old, 0, 0);
      psc_mfields_destroy(flds_base_old);
    }
    *p->flds_p = flds_base_new;
  }

  communicate_free(ctx);

  psc_balance_seed_patches(domain_old, domain_new);	//TODO required here?

  psc->mrc_domain = domain_new;
  psc_bnd_particles_check_domain(psc->bnd_particles);
  psc_bnd_check_domain(psc->bnd);
  psc_output_fields_check_bnd = true;
  mrc_domain_destroy(domain_old);
}

// FIXME, way too much duplication from the above

void
psc_balance_run(struct psc_balance *bal, struct psc *psc)
{
  if (bal->force_update == true)
  {
    bal->force_update = false;
  }
  else
  {
    if (bal->every <= 0)
      return;

    if (psc->timestep == 0 || psc->timestep % bal->every != 0)
      return;
  }

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
  double *loads = calloc(nr_patches, sizeof(*loads));
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

  free(psc->patch);
  struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_new);
  prof_stop(pr_bal_decomp_B);
  prof_start(pr_bal_decomp_C);
  //  mrc_domain_view(domain_new);
  psc_setup_patches(psc, domain_new);
  prof_stop(pr_bal_decomp_C);
  prof_start(pr_bal_decomp_D);
  free(psc_balance_comp_time_by_patch);
  psc_balance_comp_time_by_patch = calloc(nr_patches_new,// leaked the very last time
					  sizeof(*psc_balance_comp_time_by_patch));
  
  //If there are no active patches, exit here
  int n_global_patches;
  mrc_domain_get_nr_global_patches(domain_new, &n_global_patches);
  if(n_global_patches < 1) exit(0);
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
  int *nr_particles_by_patch = calloc(nr_patches, sizeof(*nr_particles_by_patch));
  for (int p = 0; p < nr_patches; p++) {
    nr_particles_by_patch[p] =
      psc_mparticles_nr_particles_by_patch(psc->particles, p);
  }
  prof_stop(pr_bal_prts_A);

  communicate_new_nr_particles(ctx, &nr_particles_by_patch);

  prof_start(pr_bal_prts_B);
  // alloc new particles
  struct psc_balance_ops *ops = psc_balance_ops(bal);

  mparticles_base_t *mprts_base_new = 
    psc_mparticles_create(mrc_domain_comm(domain_new));
  psc_mparticles_set_type(mprts_base_new, psc->prm.particles_base);
  psc_mparticles_set_domain_nr_particles(mprts_base_new, domain_new,
					 nr_particles_by_patch);
  unsigned int mp_flags;
  psc_mparticles_get_param_int(psc->particles, "flags", (int *) &mp_flags);
  psc_mparticles_set_param_int(mprts_base_new, "flags", mp_flags);

  struct psc_mparticles *mprts_old = psc_mparticles_get_as(psc->particles, ops->mprts_type, 0);
  if (mprts_old != psc->particles) { // FIXME hacky: destroy old particles early if we just got a copy
    psc_mparticles_destroy(psc->particles);
  }

  prof_start(pr_bal_prts_B1);
  psc_mparticles_setup(mprts_base_new);
  prof_stop(pr_bal_prts_B1);

  struct psc_mparticles *mprts_new = psc_mparticles_get_as(mprts_base_new, ops->mprts_type, MP_DONT_COPY);
  prof_stop(pr_bal_prts_B);
    
  // communicate particles
  psc_balance_communicate_particles(bal, ctx, mprts_old, mprts_new, nr_particles_by_patch);

  prof_start(pr_bal_prts_C);
  free(nr_particles_by_patch);

  psc_mparticles_put_as(mprts_new, mprts_base_new, 0);
  psc_mparticles_put_as(mprts_old, psc->particles, MP_DONT_COPY);

  if (mprts_old == psc->particles) {
    psc_mparticles_destroy(psc->particles);
  }
  // replace particles by redistributed ones
  psc->particles = mprts_base_new;
  prof_stop(pr_bal_prts_C);

  prof_stop(pr_bal_prts);

  // ----------------------------------------------------------------------
  // fields

  prof_start(pr_bal_flds_comm);
  prof_stop(pr_bal_flds_comm);
  /* prof_start(pr_bal_flds_A); */
  /* prof_stop(pr_bal_flds_A); */
  prof_start(pr_bal_flds_B);
  prof_stop(pr_bal_flds_B);

  prof_start(pr_bal_flds);
  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
    mfields_base_t *flds_base_old = *p->flds_p;
    // if (flds_base_old != psc->flds) {
    //   fprintf(stderr, "WARNING: not rebalancing some extra field (%p) -- expect crash!\n",
    // 	      flds_base_old);
    //   continue; // FIXME!!!
    // }
    mfields_base_t *flds_base_new;
    flds_base_new = psc_mfields_create(mrc_domain_comm(domain_new));
    psc_mfields_set_type(flds_base_new, psc_mfields_type(flds_base_old));
    psc_mfields_set_name(flds_base_new, psc_mfields_name(flds_base_old));
    psc_mfields_set_domain(flds_base_new, domain_new);
    psc_mfields_set_param_int(flds_base_new, "nr_fields", flds_base_old->nr_fields);
    psc_mfields_set_param_int(flds_base_new, "first_comp", flds_base_old->first_comp);
    psc_mfields_set_param_int3(flds_base_new, "ibn", flds_base_old->ibn);

    prof_start(pr_bal_flds_A);
    psc_mfields_setup(flds_base_new);
    prof_stop(pr_bal_flds_A);
    for (int m = flds_base_old->first_comp;
	 m < flds_base_old->first_comp + flds_base_old->nr_fields; m++) {
      const char *s = psc_mfields_comp_name(flds_base_old, m);
      if (s) {
	psc_mfields_set_comp_name(flds_base_new, m, s);
      }
    }

    prof_restart(pr_bal_flds_B);
    struct psc_mfields *flds_old =
      psc_mfields_get_as(flds_base_old, ops->mflds_type, flds_base_old->first_comp,
			 flds_base_old->first_comp + flds_base_old->nr_fields);
    struct psc_mfields *flds_new =
      psc_mfields_get_as(flds_base_new, ops->mflds_type, 0, 0);
    prof_stop(pr_bal_flds_B);

    prof_restart(pr_bal_flds_comm);
    psc_balance_communicate_fields(bal, ctx, flds_old, flds_new);
    prof_stop(pr_bal_flds_comm);

    prof_restart(pr_bal_flds_C);
    psc_mfields_put_as(flds_old, flds_base_old, 0, 0);
    psc_mfields_put_as(flds_new, flds_base_new, flds_base_new->first_comp,
		       flds_base_new->first_comp + flds_base_new->nr_fields);
    prof_stop(pr_bal_flds_C);

    psc_mfields_destroy(*p->flds_p);
    *p->flds_p = flds_base_new;
  }
  prof_stop(pr_bal_flds);
  
  communicate_free(ctx);

  psc_balance_seed_patches(domain_old, domain_new);

  // ----------------------------------------------------------------------
  // photons
  // alloc new photons
  // FIXME, will break if there are actual photons
  mphotons_t *mphotons_new = psc_mphotons_create(mrc_domain_comm(domain_new));
  psc_mphotons_set_domain(mphotons_new, domain_new);
  psc_mphotons_setup(mphotons_new);

  // replace photons by redistributed ones
  psc_mphotons_destroy(psc->mphotons);
  psc->mphotons = mphotons_new;


  psc->mrc_domain = domain_new;
  psc_bnd_particles_check_domain(psc->bnd_particles);
  psc_bnd_check_domain(psc->bnd);
  psc_output_fields_check_bnd = true;
  
  mrc_domain_destroy(domain_old);

  psc_stats_stop(st_time_balance);
}

// ----------------------------------------------------------------------
// _psc_balance_read

static void
_psc_balance_read(struct psc_balance *bal, struct mrc_io *io)
{
  int nr_patches;
  mrc_domain_get_patches(ppsc->mrc_domain, &nr_patches);
  psc_balance_comp_time_by_patch = calloc(nr_patches,// leaked the very last time
					  sizeof(*psc_balance_comp_time_by_patch));
}

// ----------------------------------------------------------------------
// _psc_balance_destroy

static void
_psc_balance_destroy(struct psc_balance *bal)
{
  free(psc_balance_comp_time_by_patch);
}

// ----------------------------------------------------------------------
// psc_balance_init

static void
psc_balance_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_balance, &psc_balance_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_balance, &psc_balance_single_ops);
}

// ======================================================================
// psc_balance class

#define VAR(x) (void *)offsetof(struct psc_balance, x)
static struct param psc_balance_descr[] = {
  { "every"            , VAR(every)            , PARAM_INT(0)            },
  { "force_update"     , VAR(force_update)     , PARAM_INT(0)            },
  { "factor_fields"    , VAR(factor_fields)    , PARAM_DOUBLE(1.)        },
  { "print_loads"      , VAR(print_loads)      , PARAM_BOOL(false)       },
  { "write_loads"      , VAR(write_loads)      , PARAM_BOOL(false)       },
  {},
};
#undef VAR

struct mrc_class_psc_balance mrc_class_psc_balance = {
  .name             = "psc_balance",
  .size             = sizeof(struct psc_balance),
  .param_descr      = psc_balance_descr,
  .init             = psc_balance_init,
  .destroy          = _psc_balance_destroy,
  .read             = _psc_balance_read,
};

