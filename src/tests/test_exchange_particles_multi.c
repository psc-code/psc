
#include "psc_testing.h"
#include "psc_bnd.h"
#include "psc_particles_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

static void
psc_test_setup_particles(struct psc *psc, int *nr_particles_by_patch, bool count_only)
{
  for (int p = 0; p < psc->nr_patches; p++) {
    struct psc_patch *patch = &ppsc->patch[p];
 
    int i = 0;
    if (p == 0) { // initially, only particles on first patch
      int *ilo = patch->off;
      int ihi[3] = { patch->off[0] + patch->ldims[0],
		     patch->off[1] + patch->ldims[1],
		     patch->off[2] + patch->ldims[2] };
      for (int iz = ilo[2]-1; iz < ihi[2]+1; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) { // xz only !!!
	  for (int ix = ilo[0]-1; ix < ihi[0]+1; ix++) {
	    if (!count_only) {
	      particles_c_t *pp = psc_mparticles_get_patch_c(psc->particles, p);
	      particle_t *p;
	      p = particles_get_one(pp, i);
	      memset(p, 0, sizeof(*p));
	      p->qni = -1.; p->mni = 1.; p->wni = 1.;
	      p->xi = (ix + .5) * ppsc->dx[0];
	      p->yi = (iy + .5) * ppsc->dx[1];
	      p->zi = (iz + .5) * ppsc->dx[2];
	      
	      p = particles_get_one(pp, i + 1);
	      memset(p, 0, sizeof(*p));
	      p->qni = 1.; p->mni = 100.; p->wni = 1.;
	      p->xi = (ix + .5) * ppsc->dx[0];
	      p->yi = (iy + .5) * ppsc->dx[1];
	      p->zi = (iz + .5) * ppsc->dx[2];
	    }
	    i += 2;
	  }
	}
      }
    }

    if (count_only) {
      nr_particles_by_patch[p] = i;
    } else {
      particles_c_t *pp = psc_mparticles_get_patch_c(psc->particles, p);
      pp->n_part = i;
    }
  }
}

// FIXME, make generic
static int
get_total_num_particles(mparticles_base_t *particles_base)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, MP_DONT_COPY);

  int nr_part = 0;
  psc_foreach_patch(ppsc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    nr_part += pp->n_part;
  }

  int total_nr_part;
  MPI_Allreduce(&nr_part, &total_nr_part, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

  psc_mparticles_put_cf(particles, particles_base);
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
  if (opt_testing_dump) {
    psc_dump_particles(psc->particles, "part-0");
  }
  int total_num_particles_before = get_total_num_particles(psc->particles);
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  if (opt_testing_dump) {
    psc_dump_particles(psc->particles, "part-1");
  }
  int total_num_particles_after = get_total_num_particles(psc->particles);
  psc_check_particles(psc->particles);
  assert(total_num_particles_before == total_num_particles_after);
  psc_destroy(psc);

  psc_testing_finalize();
}
