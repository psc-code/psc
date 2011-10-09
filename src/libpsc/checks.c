
#include "psc.h"
#include "psc_particles_as_c.h"

void
psc_check_particles(mparticles_base_t *particles_base)
{
  int fail_cnt = 0;

  mparticles_t particles;
  psc_mparticles_base_get_cf(&particles, particles_base);

  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particles_t *pp = &particles.p[p];
    f_real xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = patch->xb[d];
      xe[d] = patch->xb[d] + patch->ldims[d] * ppsc->dx[d];
    }
    
    for (int i = 0; i < pp->n_part; i++) {
      particle_t *part = particles_get_one(pp, i);
      if (part->xi < xb[0] || part->xi >= xe[0] || // FIXME xz only!
	  part->zi < xb[2] || part->zi >= xe[2]) {
	if (fail_cnt++ < 10) {
	  mprintf("FAIL: xi %g [%g:%g]\n", part->xi, xb[0], xe[0]);
	  mprintf("      zi %g [%g:%g]\n", part->zi, xb[2], xe[2]);
	}
      }
    }
  }
  assert(fail_cnt == 0);
  psc_mparticles_base_put_cf(&particles, particles_base); // FIXME, no copy-back needed
}

