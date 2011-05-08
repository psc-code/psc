
#include "psc.h"

void
psc_check_particles(mparticles_base_t *particles)
{
  int fail_cnt = 0;
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_base_t *pp = &particles->p[p];
    f_real xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = (-.5) * ppsc->dx[d] + patch->xb[d];
      xe[d] = (patch->ldims[d]-.5) * ppsc->dx[d] + patch->xb[d];
    }
    
    for (int i = 0; i < pp->n_part; i++) {
      particle_base_t *part = particles_base_get_one(pp, i);
      if (part->xi < xb[0] || part->xi > xe[0] || // FIXME xz only!
	  part->zi < xb[2] || part->zi > xe[2]) {
	if (fail_cnt++ < 10) {
	  mprintf("FAIL: xi %g [%g:%g]\n", part->xi, xb[0], xe[0]);
	  mprintf("      zi %g [%g:%g]\n", part->zi, xb[2], xe[2]);
	}
      }
    }
  }
  assert(fail_cnt == 0);
  MHERE;
}

