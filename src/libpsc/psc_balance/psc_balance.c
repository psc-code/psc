
#include "psc_balance.h"
#include "psc_case.h"
#include "psc_bnd_fields.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

LIST_HEAD(psc_mfields_base_list);

struct psc_balance {
  struct mrc_obj obj;
  int every;
  int force_update;
};

static void
psc_get_loads_initial(struct psc *psc, double *loads, int *nr_particles_by_patch)
{
  psc_foreach_patch(psc, p) {
    loads[p] = nr_particles_by_patch[p] + 1;
  }
}

static void
psc_get_loads(struct psc *psc, double *loads)
{
  mparticles_t *mparticles = psc_mparticles_get_cf(psc->particles, 0);

  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(mparticles, p);
    loads[p] = pp->n_part + 1;
  }

  psc_mparticles_put_cf(mparticles, psc->particles); // OPT, doesn't need copy back
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
    nr_patches_all_new = calloc(size, sizeof(*nr_patches_all_new));
    double loads_sum = 0.;
    for (int i = 0; i < nr_global_patches; i++) {
      loads_sum += loads_all[i];
    }
    double load_target = loads_sum / size;
    //    mprintf("target %g size %d sum %g\n", load_target, size, loads_sum);
    
    int p = 0, nr_new_patches = 0;
    double load = 0.;
    for (int i = 0; i < nr_global_patches; i++) {
      load += loads_all[i];
      nr_new_patches++;
      double next_target = (p + 1) * load_target;
      if (p < size - 1) {
	if (load > next_target) {
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
	}
      }
      // last proc takes what's left
      if (i == nr_global_patches - 1) {
	nr_patches_all_new[size - 1] = nr_new_patches;
      }
    }
    
    int pp = 0;
    for (int p = 0; p < size; p++) {
      double load = 0.;
      for (int i = 0; i < nr_patches_all_new[p]; i++) {
	load += loads_all[pp++];
	//	mprintf("  pp %d load %g : %g\n", pp-1, loads_all[pp-1], load);
      }
      //      mprintf("p %d # = %d load %g / %g : %g\n", p, nr_patches_all_new[p],
      //	      load, load_target, load - load_target);
    }
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

  if (rank == 0) {
    free(nr_patches_all);
    free(displs);
  }

  return loads_all;
}

static inline int
mpi_tag(struct mrc_patch_info *info)
{
  // This works up to 1000 patches per direction
  return info->idx3[0] + (info->idx3[1] << 10) + (info->idx3[2] << 20);
}

static void
communicate_new_nr_particles(struct mrc_domain *domain_old,
			     struct mrc_domain *domain_new, int **p_nr_particles_by_patch)
{
  MPI_Comm comm = mrc_domain_comm(domain_new);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int nr_patches_old, nr_patches_new;
  mrc_domain_get_patches(domain_old, &nr_patches_old);
  mrc_domain_get_patches(domain_new, &nr_patches_new);

  int *nr_particles_by_patch_old = *p_nr_particles_by_patch;
  int *nr_particles_by_patch_new = calloc(nr_patches_new,
					  sizeof(nr_particles_by_patch_new));

  MPI_Request *send_reqs = calloc(nr_patches_old, sizeof(*send_reqs));
  // send info from old local patches
  for (int p = 0; p < nr_patches_old; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info);
    mrc_domain_get_idx3_patch_info(domain_new, info.idx3, &info_new);
    if (info_new.rank == rank) {
      nr_particles_by_patch_new[info_new.patch] = nr_particles_by_patch_old[p];
      send_reqs[p] = MPI_REQUEST_NULL;
    } else if ( info_new.rank < 0 ) {	//This patch has been deleted
      //Issue a warning if there will be particles lost
      if(nr_particles_by_patch_old[p] > 0) printf("Warning: losing %d particles in patch deallocation\n", nr_particles_by_patch_old[p]);
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      MPI_Isend(&nr_particles_by_patch_old[p], 1, MPI_INT, info_new.rank,
		mpi_tag(&info), comm, &send_reqs[p]);
    }
  }
  // recv info for new local patches
  MPI_Request *recv_reqs = calloc(nr_patches_new, sizeof(*recv_reqs));
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if ( info_old.rank < 0 ) {
      //TODO Get number of particles
      nr_particles_by_patch_new[p] = psc_case_calc_nr_particles_in_patch(ppsc->patchmanager.currentcase, p);
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else {
      //printf("a: rank: %d tag: %d\n", info_old.rank, mpi_tag(&info));
      MPI_Irecv(&nr_particles_by_patch_new[p], 1, MPI_INT, info_old.rank, mpi_tag(&info),
		comm, &recv_reqs[p]);
    }
  }
  
  MPI_Waitall(nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
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
    particles_t *pp = psc_mparticles_get_patch(particles_new, p);
    pp->n_part = nr_particles_by_patch_new[p];
  }

  MPI_Request *send_reqs = calloc(nr_patches_old, sizeof(*send_reqs));
  // send from old local patches
  for (int p = 0; p < nr_patches_old; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info);
    mrc_domain_get_idx3_patch_info(domain_new, info.idx3, &info_new);
    if (info_new.rank == rank || info_new.rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      particles_t *pp_old = psc_mparticles_get_patch(particles_old, p);
      int nn = pp_old->n_part * (sizeof(particle_t)  / sizeof(particle_real_t));
      MPI_Isend(pp_old->particles, nn, MPI_PARTICLES_REAL, info_new.rank,
		mpi_tag(&info), comm, &send_reqs[p]);
    }
  }

  // recv for new local patches
  MPI_Request *recv_reqs = calloc(nr_patches_new, sizeof(*recv_reqs));
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (info_old.rank < 0) {
      recv_reqs[p] = MPI_REQUEST_NULL;
      //TODO Seed particles
    } else {
      particles_t *pp_new = psc_mparticles_get_patch(particles_new, p);
      int nn = pp_new->n_part * (sizeof(particle_t)  / sizeof(particle_real_t));
      MPI_Irecv(pp_new->particles, nn, MPI_PARTICLES_REAL, info_old.rank,
		mpi_tag(&info), comm, &recv_reqs[p]);
    }
  }

  // local particles
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank != rank) {
      continue;
    }

    particles_t *pp_old = psc_mparticles_get_patch(particles_old, info_old.patch);
    particles_t *pp_new = psc_mparticles_get_patch(particles_new, p);
    assert(pp_old->n_part == pp_new->n_part);
    for (int n = 0; n < pp_new->n_part; n++) {
      pp_new->particles[n] = pp_old->particles[n];
    }
  }
  
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
  
  fields_t *f_old = psc_mfields_get_patch(flds_old, 0);
  fields_t *f_new = psc_mfields_get_patch(flds_new, 0);

  MPI_Request *send_reqs = calloc(nr_patches_old, sizeof(*send_reqs));
  // send from old local patches
  for (int p = 0; p < nr_patches_old; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain_old, p, &info);
    mrc_domain_get_idx3_patch_info(domain_new, info.idx3, &info_new);
    if (info_new.rank == rank || info_new.rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      fields_t *pf_old = &f_old[p];
      int nn = fields_size(pf_old) * pf_old->nr_comp;
      int *ib = pf_old->ib;
      void *addr_old = &F3(pf_old, 0, ib[0], ib[1], ib[2]);
      MPI_Isend(addr_old, nn, MPI_FIELDS_REAL, info_new.rank,
		mpi_tag(&info), comm, &send_reqs[p]);
    }
  }

  // recv for new local patches
  MPI_Request *recv_reqs = calloc(nr_patches_new, sizeof(*recv_reqs));
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (info_old.rank < 0) {	//this patch did not exist before
      recv_reqs[p] = MPI_REQUEST_NULL;
      //Seed new data
    }
    else {
      fields_t *pf_new = &f_new[p];
      int nn = fields_size(pf_new) * pf_new->nr_comp;
      int *ib = pf_new->ib;
      void *addr_new = &F3(pf_new, 0, ib[0], ib[1], ib[2]);
      MPI_Irecv(addr_new, nn, MPI_FIELDS_REAL, info_old.rank,
		mpi_tag(&info), comm, &recv_reqs[p]);
    }
  }

  // local fields
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < nr_patches_new; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank != rank) {
      continue;
    }

    fields_t *pf_old = &f_old[info_old.patch];
    fields_t *pf_new = &f_new[p];

    assert(pf_old->nr_comp == pf_new->nr_comp);
    assert(fields_size(pf_old) == fields_size(pf_new));
    int size = fields_size(pf_old) * pf_old->nr_comp;
    int *ib = pf_new->ib;
    void *addr_new = &F3(pf_new, 0, ib[0], ib[1], ib[2]);
    void *addr_old = &F3(pf_old, 0, ib[0], ib[1], ib[2]);
    memcpy(addr_new, addr_old, size * sizeof(fields_real_t));
  }

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
    mrc_domain_get_idx3_patch_info(domain_old, info.idx3, &info_old);
    if (info_old.rank < 0)	//Patch has to be seeded
    {
      //Seed field
      double t = ppsc->timestep * ppsc->dt;
      psc_bnd_fields_setup_patch(psc_push_fields_get_bnd_fields(ppsc->push_fields), p, ppsc->flds, t);

      //Seed particles
      psc_case_init_particles_patch(ppsc->patchmanager.currentcase, p, 0);
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

  communicate_new_nr_particles(domain_old, domain_new, p_nr_particles_by_patch);

  // ----------------------------------------------------------------------
  // fields

  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
    mfields_base_t *flds_base_old = *p->flds_p;
    
    if (flds_base_old != psc->flds) {
      fprintf(stderr, "WARNING: not rebalancing some extra field -- expect crash!\n");
      continue; // FIXME!!!
    }
    mfields_base_t *flds_base_new;
    flds_base_new = psc_mfields_create(mrc_domain_comm(domain_new));
    psc_mfields_set_type(flds_base_new, psc_mfields_type(flds_base_old));
    psc_mfields_set_name(flds_base_new, "mfields");
    psc_mfields_set_domain(flds_base_new, domain_new);
    psc_mfields_set_param_int(flds_base_new, "nr_fields", flds_base_old->nr_fields);
    psc_mfields_set_param_int(flds_base_new, "first_comp", flds_base_old->first_comp);
    psc_mfields_set_param_int3(flds_base_new, "ibn", flds_base_old->ibn);
    psc_mfields_setup(flds_base_new);

    mfields_t *flds_old =
      psc_mfields_get_cf(flds_base_old, flds_base_old->first_comp,
			 flds_base_old->first_comp + flds_base_old->nr_fields);
    mfields_t *flds_new = psc_mfields_get_cf(flds_base_new, 0, 0);
    communicate_fields(domain_old, domain_new, flds_old, flds_new);
    psc_mfields_put_cf(flds_old, flds_base_old, 0, 0);
    psc_mfields_put_cf(flds_new, flds_base_new, flds_base_new->first_comp,
		       flds_base_new->first_comp + flds_base_new->nr_fields);

    psc_mfields_destroy(psc->flds);
    psc->flds = flds_base_new;
  }

  psc_balance_seed_patches(domain_old, domain_new);	//TODO required here?

  mrc_domain_destroy(domain_old);
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

  psc_stats_start(st_time_balance);
  struct mrc_domain *domain_old = psc->mrc_domain;

  int nr_patches;
  mrc_domain_get_patches(domain_old, &nr_patches);
  double *loads = calloc(nr_patches, sizeof(*loads));
  psc_get_loads(psc, loads);

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
  
  //If there are no active patches, exit here
  int n_global_patches;
  mrc_domain_get_nr_global_patches(domain_new, &n_global_patches);
  if(n_global_patches < 1) exit(0);

  int *nr_particles_by_patch = calloc(nr_patches, sizeof(*nr_particles_by_patch));
  for (int p = 0; p < nr_patches; p++) {
    nr_particles_by_patch[p] =
      psc_mparticles_nr_particles_by_patch(psc->particles, p);
  }
  communicate_new_nr_particles(domain_old, domain_new, &nr_particles_by_patch);

  // OPT: if local patches didn't change at all, no need to do anything...

  // ----------------------------------------------------------------------
  // particles

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

  psc_mparticles_put_cf(mparticles_old, psc->particles); // FIXME, don't need copy-back
  psc_mparticles_put_cf(mparticles_new, mparticles_base_new);

  // replace particles by redistributed ones
  psc_mparticles_destroy(psc->particles);
  psc->particles = mparticles_base_new;

  // ----------------------------------------------------------------------
  // fields

  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, &psc_mfields_base_list, entry, struct psc_mfields_list_entry) {
    mfields_base_t *flds_base_old = *p->flds_p;
    if (flds_base_old != psc->flds) {
      fprintf(stderr, "WARNING: not rebalancing some extra field -- expect crash!\n");
      continue; // FIXME!!!
    }
    mfields_base_t *flds_base_new;
    flds_base_new = psc_mfields_create(mrc_domain_comm(domain_new));
    psc_mfields_set_type(flds_base_new, psc->prm.fields_base);
    psc_mfields_set_name(flds_base_new, "mfields");
    psc_mfields_set_domain(flds_base_new, domain_new);
    psc_mfields_set_param_int(flds_base_new, "nr_fields", NR_FIELDS);
    psc_mfields_set_param_int3(flds_base_new, "ibn", psc->ibn);
    psc_mfields_setup(flds_base_new);

    mfields_t *flds_old = psc_mfields_get_cf(flds_base_old, 0, 12); // FIXME NR_FIELDS?
    mfields_t *flds_new = psc_mfields_get_cf(flds_base_new, 0, 0);
    communicate_fields(domain_old, domain_new, flds_old, flds_new);
    psc_mfields_put_cf(flds_old, flds_base_old, 0, 0);
    psc_mfields_put_cf(flds_new, flds_base_new, 0, 12);

    psc_mfields_destroy(psc->flds);
    psc->flds = flds_base_new;
  }
  
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


  mrc_domain_destroy(domain_old);
  psc->mrc_domain = domain_new;

  psc_stats_stop(st_time_balance);
}

// ======================================================================
// psc_balance class

#define VAR(x) (void *)offsetof(struct psc_balance, x)
static struct param psc_balance_descr[] = {
  { "every"            , VAR(every)               , PARAM_INT(0)        },
  { "force_update"     , VAR(force_update)	  , PARAM_INT(0)	},
  {},
};
#undef VAR

struct mrc_class_psc_balance mrc_class_psc_balance = {
  .name             = "psc_balance",
  .size             = sizeof(struct psc_balance),
  .param_descr      = psc_balance_descr,
};

