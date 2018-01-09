
#include "psc_testing.h"
#include "psc_bnd_particles.h"
#include "psc_particles_as_double.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

static void
psc_test_setup_particles(struct psc *psc, int *nr_particles_by_patch, bool count_only)
{
  assert(0);
#if 0
  int rank;
  MPI_Comm_rank(psc_comm(psc), &rank);

  mparticles_t mprts(psc->particles);

  for (int p = 0; p < psc->nr_patches; p++) {
    struct psc_patch *patch = &ppsc->patch[p];
 
    int i = 0;
    if (p == 0 && rank == 0) { // initially, only particles on first patch
      int *ilo = patch->off;
      int ihi[3] = { patch->off[0] + patch->ldims[0],
		     patch->off[1] + patch->ldims[1],
		     patch->off[2] + patch->ldims[2] };
      for (int iz = ilo[2]-1; iz < ihi[2]+1; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) { // xz only !!!
	  for (int ix = ilo[0]-1; ix < ihi[0]+1; ix++) {
	    i += 2;
	    if (!count_only) {
	      particle_t prt = {};
	      prt.qni = -1.; prt.mni = 1.; prt.wni = 1.;
	      prt.xi = (ix + .5) * patch->dx[0];
	      prt.yi = (iy + .5) * patch->dx[1];
	      prt.zi = (iz + .5) * patch->dx[2];
	      mprts[p].push_back(prt);
	      
	      prt.qni = 1.; prt.mni = 100.; prt.wni = 1.;
	      prt.xi = (ix + .5) * patch->dx[0];
	      prt.yi = (iy + .5) * patch->dx[1];
	      prt.zi = (iz + .5) * patch->dx[2];
	      mprts[p].push_back(prt);
	    }
	  }
	}
      }
    }

    if (count_only) {
      nr_particles_by_patch[p] = i;
    }
  }
#endif
}

// FIXME, make generic
static int
get_total_num_particles(struct psc_mparticles *particles_base)
{
  int nr_part = psc_mparticles_nr_particles(particles_base);
  int total_nr_part;
  MPI_Allreduce(&nr_part, &total_nr_part, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

  return total_nr_part;
}

// ======================================================================
// psc_test_ops_xz

static struct psc_ops psc_test_ops_xz = {
  .name             = "test",
  .size             = sizeof(struct psc_test),
  .create           = psc_test_create,
  .setup_particles  = psc_test_setup_particles,
};

int
main(int argc, char **argv)
{
  psc_testing_init(&argc, &argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops_xz);

  struct psc *psc = psc_testing_create_test_xz();
  psc_setup(psc);
  int total_num_particles_before = get_total_num_particles(psc->particles);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  int total_num_particles_after = get_total_num_particles(psc->particles);
  psc_mparticles_check(psc->particles);
  assert(total_num_particles_before == total_num_particles_after);
  psc_destroy(psc);

  psc_testing_finalize();
}
