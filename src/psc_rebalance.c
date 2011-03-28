
#include "psc.h"
#include "psc_output_fields.h"
#include "psc_bnd.h"

#include <stdlib.h>

static void
psc_get_loads(struct psc *psc, double *loads, int *nr_particles_by_patch)
{
  psc_foreach_patch(psc, p) {
    loads[p] = nr_particles_by_patch[p];
  }
}

static void
find_best_mapping(int nr_global_patches, double *loads_all,
		  int size, int *nr_patches_all_new)
{
  double loads_sum = 0.;
  for (int i = 0; i < nr_global_patches; i++) {
    loads_sum += loads_all[i];
  }
  double load_target = loads_sum / size;
    
  int p = 0, nr_new_patches = 0;
  double load = 0.;
  for (int i = 0; i < nr_global_patches; i++) {
    load += loads_all[i];
    nr_new_patches++;
    double next_target = (p + 1) * load_target;
    if (p < size - 1) {
      if (load > next_target) {
	double above_target = load - next_target;
	double below_target = next_target - (load - loads_all[i-1]);
	if (above_target > below_target) {
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
    }
    printf("p %d # = %d load %g / %g : %g\n", p, nr_patches_all_new[p],
	   load, load_target, load - load_target);
  }
}

void
psc_rebalance_run(struct psc *psc, int **p_nr_particles_by_patch)
{
  struct mrc_domain *domain = psc->mrc_domain;
  int *nr_particles_by_patch = *p_nr_particles_by_patch;

  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);
  double *loads = calloc(nr_patches, sizeof(*loads));
  psc_get_loads(psc, loads, nr_particles_by_patch);

  MPI_Comm comm = mrc_domain_comm(domain);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // gather nr_patches for all procs on proc 0
  int *nr_patches_all = NULL, *displs = NULL;
  int nr_global_patches = -1;
  double *loads_all = NULL;
  if (rank == 0) {
    nr_patches_all = calloc(size, sizeof(*nr_patches_all));
    displs = calloc(size, sizeof(*displs));
  }
  MPI_Gather(&nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);
  // gather loads for all patches on proc 0
  if (rank == 0) {
    nr_global_patches = 0;
    for (int i = 0; i < size; i++) {
      displs[i] = nr_global_patches;
      nr_global_patches += nr_patches_all[i];
    }
    loads_all = calloc(nr_global_patches, sizeof(*loads_all));
  }
  MPI_Gatherv(loads, nr_patches, MPI_DOUBLE, loads_all, nr_patches_all, displs,
	      MPI_DOUBLE, 0, comm);

  int *nr_patches_all_new = calloc(size, sizeof(*nr_patches_all_new));
  if (rank == 0) {
    find_best_mapping(nr_global_patches, loads_all, size, nr_patches_all_new);
  }
  MPI_Bcast(nr_patches_all_new, size, MPI_INT, 0, comm);

  if (rank == 0) {
    free(loads_all);
    free(nr_patches_all);
    free(displs);
  }
  
  free(loads);

  free(psc->patch);
  struct mrc_domain *domain_new = psc_setup_mrc_domain(psc, nr_patches_all_new[rank]);

  free(nr_patches_all_new);

  mrc_domain_view(domain_new);

  int *nr_particles_by_patch_new = calloc(psc->nr_patches,
					  sizeof(nr_particles_by_patch_new));

  MPI_Request *send_reqs = calloc(nr_patches, sizeof(*send_reqs));
  // send info from old local patches
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch_info info, info_new;
    mrc_domain_get_local_patch_info(domain, p, &info);
    mrc_domain_get_global_patch_info(domain_new, info.global_patch, &info_new);
    if (info_new.rank == rank) {
      nr_particles_by_patch_new[info_new.patch] = nr_particles_by_patch[p];
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      MPI_Isend(&nr_particles_by_patch[p], 1, MPI_INT, info_new.rank, info.global_patch,
		comm, &send_reqs[p]);
    }
  }
  // recv info for new local patches
  MPI_Request *recv_reqs = calloc(psc->nr_patches, sizeof(*recv_reqs));
  for (int p = 0; p < psc->nr_patches; p++) {
    struct mrc_patch_info info, info_old;
    mrc_domain_get_local_patch_info(domain_new, p, &info);
    mrc_domain_get_global_patch_info(domain, info.global_patch, &info_old);
    if (info_old.rank == rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else {
      MPI_Irecv(&nr_particles_by_patch_new[p], 1, MPI_INT, info_old.rank, info.global_patch,
		comm, &recv_reqs[p]);
    }
  }
  
  MPI_Waitall(nr_patches, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(psc->nr_patches, recv_reqs, MPI_STATUSES_IGNORE);
  free(send_reqs);
  free(recv_reqs);
  
  mrc_domain_destroy(psc->mrc_domain);
  psc->mrc_domain = domain_new;

  psc_output_fields_destroy(psc->output_fields);
  psc->output_fields = psc_output_fields_create(MPI_COMM_WORLD);
  psc_output_fields_set_from_options(psc->output_fields);
  psc_output_fields_setup(psc->output_fields);

  psc_bnd_destroy(psc->bnd);
  psc->bnd = psc_bnd_create(MPI_COMM_WORLD);
  psc_bnd_set_from_options(psc->bnd);
  psc_bnd_setup(psc->bnd);

  free(*p_nr_particles_by_patch);
  *p_nr_particles_by_patch = nr_particles_by_patch_new;
}

