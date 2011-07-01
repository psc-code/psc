
#include "psc_generic_c.h"
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>

// test w/ C data structure, single precision 

static void
do_genc_push_part_yz_a(particles_t *pp)
{
  creal dt = ppsc->dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;

  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);

    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}

void
psc_push_particles_generic_c_push_yz_a(struct psc_push_particles *push,
				       mparticles_base_t *particles_base,
				       mfields_base_t *flds_base)
{
  mparticles_t particles;
  mparticles_get(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz_a", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    do_genc_push_part_yz_a(&particles.p[p]);
  }
  prof_stop(pr);

  mparticles_put(&particles, particles_base);
}
