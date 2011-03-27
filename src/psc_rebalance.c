
#include "psc.h"

#include <stdlib.h>

static void
psc_get_loads(struct psc *psc, double *loads)
{
  mparticles_base_t *mparticles = &psc->particles;

  psc_foreach_patch(psc, p) {
    loads[p] = mparticles->p[p].n_part;
  }
}

static void
find_best_assignment(int nr_global_patches, double *loads_all,
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
    } else { // last proc takes what's left
      if (i == nr_global_patches - 1) {
	nr_patches_all_new[p] = nr_new_patches;
      }
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
psc_rebalance_run(struct psc *psc)
{
  struct mrc_domain *domain = psc->mrc_domain;

  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);
  double *loads = calloc(nr_patches, sizeof(*loads));
  psc_get_loads(psc, loads);

  MPI_Comm comm = mrc_domain_comm(domain);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int *nr_patches_all = NULL, *displs = NULL;
  int nr_global_patches = -1;
  double *loads_all = NULL;
  if (rank == 0) {
    nr_patches_all = calloc(size, sizeof(*nr_patches_all));
    displs = calloc(size, sizeof(*displs));
  }
  MPI_Gather(&nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);
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

  int *nr_patches_all_new;
  if (rank == 0) {
    nr_patches_all_new = calloc(size, sizeof(*nr_patches_all_new));
    find_best_assignment(nr_global_patches, loads_all, size, nr_patches_all_new);
  }

  if (rank == 0) {
    free(nr_patches_all_new);
    free(loads_all);
    free(nr_patches_all);
    free(displs);
  }
  
  free(loads);

  MPI_Barrier(MPI_COMM_WORLD);
  //  assert(0);
}

