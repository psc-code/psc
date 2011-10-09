
#include "psc_testing.h"
#include "psc_bnd.h"
#include "psc_particles_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

void
setup_particles(mparticles_base_t *particles_base)
{
  mparticles_t *particles = psc_mparticles_base_get_cf(particles_base);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // FIXME, realloc
  // only particles on proc 0 / patch 0, but some are out of bounds.
  // particles are right on nodes of the grid, but as far as
  // particle -> cell mapping goes, that's right in the center of the
  // particle ordering cell, for reasons that should be laid out more clearly.
  // for the old ordering, nodes aren't good because they're indeterminate
  // (could go either way), so let's shift them a bit so we get a unique answer
  // we can check.

  psc_foreach_patch(ppsc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    if (rank != 0 || p != 0) {
      pp->n_part = 0;
      continue;
    }

    struct psc_patch *patch = &ppsc->patch[p];
    int *ilo = patch->off;
    int ihi[3] = { patch->off[0] + patch->ldims[0],
		   patch->off[1] + patch->ldims[1],
		   patch->off[2] + patch->ldims[2] };
    int i = 0;
    
    for (int iz = ilo[2]-1; iz < ihi[2]+1; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) { // xz only !!!
	for (int ix = ilo[0]-1; ix < ihi[0]+1; ix++) {
	  particle_t *p;
	  p = particles_get_one(pp, i++);
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix + .51) * ppsc->dx[0];
	  p->yi = (iy + .51) * ppsc->dx[1];
	  p->zi = (iz + .51) * ppsc->dx[2];
	  
	  p = particles_get_one(pp, i++);
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix + .49) * ppsc->dx[0];
	  p->yi = (iy + .49) * ppsc->dx[1];
	  p->zi = (iz + .49) * ppsc->dx[2];
	}
      }
    }
    pp->n_part = i;
  }
  psc_mparticles_base_put_cf(particles, particles_base);
}

// FIXME, make generic
static int
get_total_num_particles(mparticles_base_t *particles_base)
{
  mparticles_t *particles = psc_mparticles_base_get_cf(particles_base);

  int nr_part = 0;
  psc_foreach_patch(ppsc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    nr_part += pp->n_part;
  }

  int total_nr_part;
  MPI_Allreduce(&nr_part, &total_nr_part, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

  psc_mparticles_base_put_cf(particles, particles_base);
  return total_nr_part;
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  // test psc_exchange_particles()

  struct psc_case *_case = psc_create_test_xz();
  psc_bnd_set_type(ppsc->bnd, "c");
  psc_case_setup(_case);
  mparticles_base_t *particles = ppsc->particles;
  setup_particles(particles);
  //  psc_dump_particles("part-0");
  int total_num_particles_before = get_total_num_particles(particles);
  psc_bnd_exchange_particles(ppsc->bnd, particles);
  //  psc_dump_particles("part-1");
  int total_num_particles_after = get_total_num_particles(particles);
  psc_check_particles(particles);
  assert(total_num_particles_before == total_num_particles_after);
  psc_case_destroy(_case);

  prof_print();

  MPI_Finalize();
}
