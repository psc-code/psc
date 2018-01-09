
#include "psc_testing.h"
#include "psc_bnd.h"
#include "psc_particles_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#if 0
void
setup_particles(struct psc_mparticles *particles_base)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, 0);
  // FIXME, realloc
  // only particles on proc 1, but some are out of bounds.
  // particles are right on nodes of the grid, but as far as
  // particle -> cell mapping goes, that's right in the center of the
  // particle ordering cell, for reasons that should be laid out more clearly.
  // for the old ordering, nodes aren't good because they're indeterminate
  // (could go either way), so let's shift them a bit so we get a unique answer
  // we can check.
  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);
  if (rank == 0) {
    struct psc_patch *patch = &ppsc->patch[0];
    int *ilo = patch->off;
    int ihi[3] = { patch->off[0] + patch->ldims[0],
		   patch->off[1] + patch->ldims[1],
		   patch->off[2] + patch->ldims[2] };
    int i = 0;

    for (int iz = ilo[2]-1; iz < ihi[2]+1; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) { // xz only !!!
	for (int ix = ilo[0]-1; ix < ihi[0]+1; ix++) {
	  particle_t *p = prts[i++];
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix + .51) * ppsc->dx[0];
	  p->yi = (iy + .51) * ppsc->dx[1];
	  p->zi = (iz + .51) * ppsc->dx[2];

	  p = prts[i++];
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix + .49) * ppsc->dx[0];
	  p->yi = (iy + .49) * ppsc->dx[1];
	  p->zi = (iz + .49) * ppsc->dx[2];
	}
      }
    }
    prts->n_part = i;
  } else {
    prts->n_part = 0;
  }
  psc_particles_put_as(prts, prts_base, 0);
}

static void
check_particles_old_xz(struct psc_mparticles *particles_base)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, 0);

  struct psc_patch *patch = &ppsc->patch[0];
  int *ilo = patch->off;
  int ihi[3] = { patch->off[0] + patch->ldims[0],
		 patch->off[1] + patch->ldims[1],
		 patch->off[2] + patch->ldims[2] };
  f_real xb[3], xe[3];
  
  // These boundaries are, I believe, what the Fortran code guarantees.
  // The rest of the code works under these conditions, but I don't like
  // not having one strict owner for each location of a particle in the domain
  // might be.
  // OTOH, there might be a point to have valid regions overlap so that particles
  // gyrating near a boundary don't keep getting bounced between processors --
  // but that's the opposite of what's happening here.
  for (int d = 0; d < 3; d++) {
    xb[d] = (ilo[d]-2) * ppsc->dx[d];
    xe[d] = (ihi[d]+1) * ppsc->dx[d];
  }

  int fail_cnt = 0;
  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);
  for (int i = 0; i < prts->n_part; i++) {
    particle_t *p = prts[i];
    if (p->xi < xb[0] || p->xi > xe[0] ||
	p->zi < xb[2] || p->zi > xe[2]) {
      if (fail_cnt++ < 10) {
	printf("FAIL old: xi %g [%g:%g]\n", p->xi, xb[0], xe[0]);
	printf("          zi %g [%g:%g]\n", p->zi, xb[2], xe[2]);
      }
    }
  }
  assert(fail_cnt == 0);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

static void
check_particles(struct psc_mparticles *particles_base)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, 0);

  struct psc_patch *patch = &ppsc->patch[0];
  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);
  int *ilo = patch->off;
  int ihi[3] = { patch->off[0] + patch->ldims[0],
		 patch->off[1] + patch->ldims[1],
		 patch->off[2] + patch->ldims[2] };
  f_real xb[3], xe[3];

  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.
  
  for (int d = 0; d < 3; d++) {
    xb[d] = ilo[d] * ppsc->dx[d];
    xe[d] = ihi[d] * ppsc->dx[d];
  }

  int fail_cnt = 0;
  for (int i = 0; i < prts->n_part; i++) {
    particle_t *p = prts[i];
    if (p->xi < xb[0] || p->xi > xe[0] ||
	p->zi < xb[2] || p->zi > xe[2]) {
      if (fail_cnt++ < 10) {
	printf("FAIL: xi %g [%g:%g]\n", p->xi, xb[0], xe[0]);
	printf("      zi %g [%g:%g]\n", p->zi, xb[2], xe[2]);
      }
    }
  }
  assert(fail_cnt == 0);
  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

static int
get_total_num_particles(struct psc_mparticles *particles_base)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, 0);

  int total_num_part;

  MPI_Allreduce(&prts_base->n_part, &total_num_part, 1, MPI_INT, MPI_SUM,
		psc_mparticles_comm(particles_base));

  return total_num_part;
}
#endif

int
main(int argc, char **argv)
{
#if 0
  psc_testing_init(&argc, &argv);

  // test psc_exchange_particles()

  struct psc_case *_case = psc_create_test_xz();
  psc_bnd_set_type(ppsc->bnd, "fortran");
  psc_case_setup(_case);

  struct psc_mparticles *particles = ppsc->particles;
  setup_particles(particles);
  int total_num_particles_before = get_total_num_particles(particles);
  psc_bnd_exchange_particles(ppsc->bnd, particles);
  int total_num_particles_after = get_total_num_particles(particles);
  check_particles_old_xz(particles);
  assert(total_num_particles_before == total_num_particles_after);
  psc_case_destroy(_case);

  _case = psc_create_test_xz();
  psc_bnd_set_type(ppsc->bnd, "c");
  psc_case_setup(_case);
  particles = ppsc->particles;

  setup_particles(particles);
  total_num_particles_before = get_total_num_particles(particles);
  psc_bnd_exchange_particles(ppsc->bnd, particles);
  total_num_particles_after = get_total_num_particles(particles);
  check_particles(particles);
  assert(total_num_particles_before == total_num_particles_after);
  psc_case_destroy(_case);

  psc_testing_finalize();
#endif
}
