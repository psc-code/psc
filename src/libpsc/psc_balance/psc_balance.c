
#include "psc_balance.h"
#include "psc_bnd_fields.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_particles_as_double.h" // FIXME, hardcoded
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG_BALANCE

LIST_HEAD(psc_mfields_base_list);

double *psc_balance_comp_time_by_patch;

struct psc_balance {
  struct mrc_obj obj;
  int every;
  int force_update;
  double factor_fields;
};

static double
capability_default(int p)
{
  return 1.;
}

static double __unused
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
    }
  }
}

static int
find_best_mapping(struct mrc_domain *domain, int nr_global_patches, double *loads_all)
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
#ifdef DEBUG_BALANCE
    mprintf("loads_sum %g capability_sum %g load_target %g\n",
	    loads_sum, capability_sum, load_target);
#endif

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
    
#ifdef DEBUG_BALANCE
    int pp = 0;
    double min_diff = 0, max_diff = 0;
    for (int p = 0; p < size; p++) {
      double load = 0.;
      for (int i = 0; i < nr_patches_all_new[p]; i++) {
	load += loads_all[pp++];
	//mprintf("  pp %d load %g : %g\n", pp-1, loads_all[pp-1], load);
      }
      double diff = load - load_target * capability[p];
      mprintf("p %d # = %d load %g / %g : diff %g\n", p, nr_patches_all_new[p],
	      load, load_target * capability[p], diff);
      if (diff < min_diff) {
	min_diff = diff;
      }
      if (diff > max_diff) {
	max_diff = diff;
      }
    }
    mprintf("achieved target %g (%g %% -- %g %%)\n", load_target,
	    100 * min_diff / load_target, 100 * max_diff / load_target);
	      
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
#endif
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
	  if(ppsc->use_dynamic_patches)
	  {
		  int n_old_patches = *p_nr_global_patches;
		  *p_nr_global_patches = bitfield3d_count_bits_set(ppsc->patchmanager.activepatches);
		  
		  loads_all = calloc(MAX(n_old_patches, *p_nr_global_patches), sizeof(*loads_all));
		  
		  for(int i=n_old_patches; i<*p_nr_global_patches; ++i)
		  {
			  loads_all[i] = 1.;	//TODO Better assumption? Like take the median or sth alike...
		  }
	  }
	  else
	  {
		  loads_all = calloc(*p_nr_global_patches, sizeof(*loads_all));
	  }
  }
  MPI_Gatherv(loads, nr_patches, MPI_DOUBLE, loads_all, nr_patches_all, displs,
	      MPI_DOUBLE, 0, comm);

#ifdef DEBUG_BALANCE
  if (rank == 0) {
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
    free(nr_patches_all);
    free(displs);
  }
#endif

  return loads_all;
}

struct send_info {
  int rank;
  int patch;
};

struct recv_info {
  int rank;
  int patch;
};

struct communicate_ctx {
  MPI_Comm comm;
  int mpi_rank;
  int mpi_size;
  int nr_patches_old;
  int nr_patches_new;
  struct send_info *send_info; // by old patch on this proc
  struct recv_info *recv_info; // by new patch on this proc
};

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
}

static void
communicate_free(struct communicate_ctx *ctx)
{
  free(ctx->send_info);
  free(ctx->recv_info);
}

static void
communicate_new_nr_particles(struct communicate_ctx *ctx, int **p_nr_particles_by_patch)
{
  // send info from old local patches

  int *nr_particles_by_patch_old = *p_nr_particles_by_patch;
  int *nr_particles_by_patch_new = calloc(ctx->nr_patches_new,
					  sizeof(nr_particles_by_patch_new));
  MPI_Request *send_reqs = calloc(ctx->nr_patches_old, sizeof(*send_reqs));
  int *nr_patches_new_by_rank = calloc(ctx->mpi_size, sizeof(*nr_patches_new_by_rank));
  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int new_rank = ctx->send_info[p].rank;
    if (new_rank == ctx->mpi_rank) {
      nr_particles_by_patch_new[ctx->send_info[p].patch] = nr_particles_by_patch_old[p];
      send_reqs[p] = MPI_REQUEST_NULL;
    } else if (new_rank < 0 ) {	//This patch has been deleted
      //Issue a warning if there will be particles lost
      if(nr_particles_by_patch_old[p] > 0) printf("Warning: losing %d particles in patch deallocation\n", nr_particles_by_patch_old[p]);
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      int tag = nr_patches_new_by_rank[new_rank]++;
      //      mprintf("send -> %d (tags %d %d) %d\n", info_new.rank, mpi_tag(&info), tag, info_new.global_patch);
      MPI_Isend(&nr_particles_by_patch_old[p], 1, MPI_INT, new_rank,
		tag, ctx->comm, &send_reqs[p]);
    }
  }
  free(nr_patches_new_by_rank);

  // recv info for new local patches

  int *nr_patches_old_by_rank = calloc(ctx->mpi_size, sizeof(*nr_patches_old_by_rank));
  MPI_Request *recv_reqs = calloc(ctx->nr_patches_new, sizeof(*recv_reqs));
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int old_rank = ctx->recv_info[p].rank;
    if (old_rank == ctx->mpi_rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (old_rank < 0) {
      //TODO Get number of particles
      assert(0);
      //nr_particles_by_patch_new[p] = psc_case_calc_nr_particles_in_patch(ppsc->patchmanager.currentcase, p);
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else {
      //printf("a: rank: %d tag: %d\n", old_rank, mpi_tag(&info));
      int tag = nr_patches_old_by_rank[old_rank]++;
      //      mprintf("recv <- %d (tags %d %d)\n", info_old.rank, mpi_tag(&info), tag);
      MPI_Irecv(&nr_particles_by_patch_new[p], 1, MPI_INT, old_rank,
		tag, ctx->comm, &recv_reqs[p]);
    }
  }
  free(nr_patches_old_by_rank);

  MPI_Waitall(ctx->nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ctx->nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
  free(send_reqs);
  free(recv_reqs);

  free(*p_nr_particles_by_patch);
  *p_nr_particles_by_patch = nr_particles_by_patch_new;
}

static void
communicate_particles(struct mrc_domain *domain_old, struct mrc_domain *domain_new,
		      mparticles_t *particles_old, mparticles_t *particles_new,
		      int *nr_particles_by_patch_new)
{
  MPI_Comm comm = mrc_domain_comm(domain_new);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int nr_patches_old, nr_patches_new;
  mrc_domain_get_patches(domain_old, &nr_patches_old);
  mrc_domain_get_patches(domain_new, &nr_patches_new);
  
  for (int p = 0; p < nr_patches_new; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles_new, p);
    prts->n_part = nr_particles_by_patch_new[p];
  }

  MPI_Request *send_reqs = calloc(nr_patches_old, sizeof(*send_reqs));
  int *nr_patches_new_by_rank = calloc(size, sizeof(*nr_patches_new_by_rank));
  // send from old local patches
  for (int p = 0; p < nr_patches_old; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_new, info.level, info.idx3, &info_new);
    if (info_new.rank == rank || info_new.rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      struct psc_particles *pp_old = psc_mparticles_get_patch(particles_old, p);
      struct psc_particles_double *c_old = psc_particles_double(pp_old);
      int nn = pp_old->n_part * (sizeof(particle_t)  / sizeof(particle_real_t));
      int tag = nr_patches_new_by_rank[info_new.rank];
      MPI_Isend(c_old->particles, nn, MPI_PARTICLES_REAL, info_new.rank,
		tag, comm, &send_reqs[p]);
    }
  }
  free(nr_patches_new_by_rank);

  // recv for new local patches
  MPI_Request *recv_reqs = calloc(nr_patches_new, sizeof(*recv_reqs));
  int *nr_patches_old_by_rank = calloc(size, sizeof(*nr_patches_new_by_rank));
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_old, info.level, info.idx3, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (info_old.rank < 0) {
      recv_reqs[p] = MPI_REQUEST_NULL;
      //TODO Seed particles
    } else {
      struct psc_particles *pp_new = psc_mparticles_get_patch(particles_new, p);
      struct psc_particles_double *c_new = psc_particles_double(pp_new);
      int nn = pp_new->n_part * (sizeof(particle_t)  / sizeof(particle_real_t));
      int tag = nr_patches_old_by_rank[info_old.rank];
      MPI_Irecv(c_new->particles, nn, MPI_PARTICLES_REAL, info_old.rank,
		tag, comm, &recv_reqs[p]);
    }
  }
  free(nr_patches_old_by_rank);

  static int pr;
  if (!pr) {
    pr = prof_register("bal prts local", 1., 0, 0);
  }

  prof_start(pr);
  // local particles
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_old, info.level, info.idx3, &info_old);
    if (info_old.rank != rank) {
      continue;
    }

    struct psc_particles *pp_old = psc_mparticles_get_patch(particles_old, info_old.patch);
    struct psc_particles_double *c_old = psc_particles_double(pp_old);
    struct psc_particles *pp_new = psc_mparticles_get_patch(particles_new, p);
    struct psc_particles_double *c_new = psc_particles_double(pp_new);
    assert(pp_old->n_part == pp_new->n_part);
#if 0
    for (int n = 0; n < pp_new->n_part; n++) {
      c_new->particles[n] = c_old->particles[n];
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
  prof_stop(pr);
  
  MPI_Waitall(nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
  free(send_reqs);
  free(recv_reqs);
}

static void
communicate_fields(struct mrc_domain *domain_old, struct mrc_domain *domain_new,
		   mfields_t *flds_old, mfields_t *flds_new)
{
  MPI_Comm comm = mrc_domain_comm(domain_new);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int nr_patches_old, nr_patches_new;
  mrc_domain_get_patches(domain_old, &nr_patches_old);
  mrc_domain_get_patches(domain_new, &nr_patches_new);
	
	//HACK: Don't communicate output fields if they don't correspond to the domain
	//This is needed e.g. for the boosted output which handles its MPI communication internally
	//printf("Field: %s\n", flds->f[0].name);
	
	if(nr_patches_old != flds_old->nr_patches /* || strncmp(flds->f[0].name, "lab", 3) == 0 */) return;
	
	
  assert(nr_patches_old == flds_old->nr_patches);
  assert(nr_patches_old > 0);
  
  // send from old local patches
  MPI_Request *send_reqs = calloc(nr_patches_old, sizeof(*send_reqs));
  int *nr_patches_new_by_rank = calloc(size, sizeof(*nr_patches_new_by_rank));
  for (int p = 0; p < nr_patches_old; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_new, info.level, info.idx3, &info_new);
    if (info_new.rank == rank || info_new.rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      fields_t *pf_old = psc_mfields_get_patch(flds_old, p);
      int nn = psc_fields_size(pf_old) * pf_old->nr_comp;
      int *ib = pf_old->ib;
      void *addr_old = &F3(pf_old, 0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_new_by_rank[info_new.rank];
      MPI_Isend(addr_old, nn, MPI_FIELDS_REAL, info_new.rank,
		tag, comm, &send_reqs[p]);
    }
  }
  free(nr_patches_new_by_rank);

  // recv for new local patches
  MPI_Request *recv_reqs = calloc(nr_patches_new, sizeof(*recv_reqs));
  int *nr_patches_old_by_rank = calloc(size, sizeof(*nr_patches_new_by_rank));
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_old, info.level, info.idx3, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (info_old.rank < 0) {	//this patch did not exist before
      recv_reqs[p] = MPI_REQUEST_NULL;
      //Seed new data
    }
    else {
      fields_t *pf_new = psc_mfields_get_patch(flds_new, p);
      int nn = psc_fields_size(pf_new) * pf_new->nr_comp;
      int *ib = pf_new->ib;
      void *addr_new = &F3(pf_new, 0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_old_by_rank[info_old.rank];
      MPI_Irecv(addr_new, nn, MPI_FIELDS_REAL, info_old.rank,
		tag, comm, &recv_reqs[p]);
    }
  }
  free(nr_patches_old_by_rank);

  static int pr;
  if (!pr) {
    pr = prof_register("bal flds local", 1., 0, 0);
  }

  prof_start(pr);
  // local fields
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_level_idx3_patch_info(domain_old, info.level, info.idx3, &info_old);
    if (info_old.rank != rank) {
      continue;
    }

    fields_t *pf_old = psc_mfields_get_patch(flds_old, info_old.patch);
    fields_t *pf_new = psc_mfields_get_patch(flds_new, p);

    assert(pf_old->nr_comp == pf_new->nr_comp);
    assert(psc_fields_size(pf_old) == psc_fields_size(pf_new));
    int size = psc_fields_size(pf_old) * pf_old->nr_comp;
    int *ib = pf_new->ib;
    void *addr_new = &F3(pf_new, 0, ib[0], ib[1], ib[2]);
    void *addr_old = &F3(pf_old, 0, ib[0], ib[1], ib[2]);
    memcpy(addr_new, addr_old, size * sizeof(fields_real_t));
  }
  prof_stop(pr);

  MPI_Waitall(nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
  free(send_reqs);
  free(recv_reqs);
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

void
psc_balance_initial(struct psc_balance *bal, struct psc *psc,
		    int **p_nr_particles_by_patch)
{
  /*if (bal->every <= 0)
    return;*/

  struct mrc_domain *domain_old = psc->mrc_domain;

  int nr_patches;
  mrc_domain_get_patches(domain_old, &nr_patches);
  double *loads = calloc(nr_patches, sizeof(*loads));
  psc_get_loads_initial(psc, loads, *p_nr_particles_by_patch);

  int nr_global_patches;
  double *loads_all = gather_loads(domain_old, loads, nr_patches,
				   &nr_global_patches);
  free(loads);

  int nr_patches_new = find_best_mapping(domain_old, nr_global_patches,
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

  communicate_free(ctx);

  // ----------------------------------------------------------------------
  // fields

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

    mfields_t *flds_old =
      psc_mfields_get_cf(flds_base_old, flds_base_old->first_comp,
			 flds_base_old->first_comp + flds_base_old->nr_fields);
    mfields_t *flds_new = psc_mfields_get_cf(flds_base_new, 0, 0);
    communicate_fields(domain_old, domain_new, flds_old, flds_new);
    psc_mfields_put_cf(flds_old, flds_base_old, 0, 0);
    psc_mfields_put_cf(flds_new, flds_base_new, flds_base_new->first_comp,
		       flds_base_new->first_comp + flds_base_new->nr_fields);

    psc_mfields_destroy(*p->flds_p);
    *p->flds_p = flds_base_new;
  }

  psc_balance_seed_patches(domain_old, domain_new);	//TODO required here?

  // FIXME mrc_domain_destroy(domain_old);
  psc->mrc_domain = domain_new;
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
    pr_bal_prts, pr_bal_flds, pr_bal_flds_comm;
  if (!pr_bal_gather) {
    pr_bal_gather = prof_register("bal gather", 1., 0, 0);
    pr_bal_decomp_A = prof_register("bal decomp A", 1., 0, 0);
    pr_bal_decomp_B = prof_register("bal decomp B", 1., 0, 0);
    pr_bal_decomp_C = prof_register("bal decomp C", 1., 0, 0);
    pr_bal_decomp_D = prof_register("bal decomp D", 1., 0, 0);
    pr_bal_prts = prof_register("bal prts", 1., 0, 0);
    pr_bal_flds = prof_register("bal flds", 1., 0, 0);
    pr_bal_flds_comm = prof_register("bal flds comm", 1., 0, 0);
  }

  MPI_Barrier(psc_comm(psc));
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

  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_decomp_A);
  int nr_patches_new = find_best_mapping(domain_old, nr_global_patches,
					 loads_all);
  prof_stop(pr_bal_decomp_A);

  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_decomp_B);
  free(loads_all);

  free(psc->patch);
  struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_new);
  prof_stop(pr_bal_decomp_B);
  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_decomp_C);
  //  mrc_domain_view(domain_new);
  psc_setup_patches(psc, domain_new);
  prof_stop(pr_bal_decomp_C);
  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_decomp_D);
  free(psc_balance_comp_time_by_patch);
  psc_balance_comp_time_by_patch = calloc(nr_patches_new,
					  sizeof(*psc_balance_comp_time_by_patch));
  
  //If there are no active patches, exit here
  int n_global_patches;
  mrc_domain_get_nr_global_patches(domain_new, &n_global_patches);
  if(n_global_patches < 1) exit(0);
  prof_stop(pr_bal_decomp_D);

  // OPT: if local patches didn't change at all, no need to do anything...

  // ----------------------------------------------------------------------
  // particles

  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_prts);
  int *nr_particles_by_patch = calloc(nr_patches, sizeof(*nr_particles_by_patch));
  for (int p = 0; p < nr_patches; p++) {
    nr_particles_by_patch[p] =
      psc_mparticles_nr_particles_by_patch(psc->particles, p);
  }

  struct communicate_ctx _ctx, *ctx = &_ctx;
  communicate_setup(ctx, domain_old, domain_new);

  communicate_new_nr_particles(ctx, &nr_particles_by_patch);

  communicate_free(ctx);

  // alloc new particles
  mparticles_base_t *mparticles_base_new = 
    psc_mparticles_create(mrc_domain_comm(domain_new));
  psc_mparticles_set_type(mparticles_base_new, psc->prm.particles_base);
  psc_mparticles_set_domain_nr_particles(mparticles_base_new, domain_new,
					      nr_particles_by_patch);
  unsigned int mp_flags;
  psc_mparticles_get_param_int(psc->particles, "flags", (int *) &mp_flags);
  psc_mparticles_set_param_int(mparticles_base_new, "flags", mp_flags);
  psc_mparticles_setup(mparticles_base_new);

  mparticles_t *mparticles_new = psc_mparticles_get_cf(mparticles_base_new, MP_DONT_COPY);
  mparticles_t *mparticles_old = psc_mparticles_get_cf(psc->particles, 0);
    
  // communicate particles
  communicate_particles(domain_old, domain_new, 
			mparticles_old, mparticles_new, nr_particles_by_patch);
  free(nr_particles_by_patch);

  psc_mparticles_put_cf(mparticles_old, psc->particles, MP_DONT_COPY);
  psc_mparticles_put_cf(mparticles_new, mparticles_base_new, 0);

  // replace particles by redistributed ones
  psc_mparticles_destroy(psc->particles);
  psc->particles = mparticles_base_new;

  prof_stop(pr_bal_prts);

  // ----------------------------------------------------------------------
  // fields

  MPI_Barrier(psc_comm(psc));
  prof_start(pr_bal_flds_comm);
  prof_stop(pr_bal_flds_comm);
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
    psc_mfields_setup(flds_base_new);
    for (int m = flds_base_old->first_comp;
	 m < flds_base_old->first_comp + flds_base_old->nr_fields; m++) {
      const char *s = psc_mfields_comp_name(flds_base_old, m);
      if (s) {
	psc_mfields_set_comp_name(flds_base_new, m, s);
      }
    }

    mfields_t *flds_old =
      psc_mfields_get_cf(flds_base_old, flds_base_old->first_comp,
			 flds_base_old->first_comp + flds_base_old->nr_fields);
    mfields_t *flds_new = psc_mfields_get_cf(flds_base_new, 0, 0);
    prof_restart(pr_bal_flds_comm);
    communicate_fields(domain_old, domain_new, flds_old, flds_new);
    prof_stop(pr_bal_flds_comm);
    psc_mfields_put_cf(flds_old, flds_base_old, 0, 0);
    psc_mfields_put_cf(flds_new, flds_base_new, flds_base_new->first_comp,
		       flds_base_new->first_comp + flds_base_new->nr_fields);

    psc_mfields_destroy(*p->flds_p);
    *p->flds_p = flds_base_new;
  }
  prof_stop(pr_bal_flds);
  
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


  // FIXME mrc_domain_destroy(domain_old);
  psc->mrc_domain = domain_new;

  psc_stats_stop(st_time_balance);
}

// ======================================================================
// psc_balance class

#define VAR(x) (void *)offsetof(struct psc_balance, x)
static struct param psc_balance_descr[] = {
  { "every"            , VAR(every)               , PARAM_INT(0)        },
  { "force_update"     , VAR(force_update)	  , PARAM_INT(0)	},
  { "factor_fields"    , VAR(factor_fields)       , PARAM_DOUBLE(1.)    },
  {},
};
#undef VAR

struct mrc_class_psc_balance mrc_class_psc_balance = {
  .name             = "psc_balance",
  .size             = sizeof(struct psc_balance),
  .param_descr      = psc_balance_descr,
};

