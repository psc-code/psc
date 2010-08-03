
#include "psc_generic_c.h"

#include "util/profile.h"
#include <math.h>
#include <stdlib.h>

// test w/ C data structure, single precision 

static void
do_genc_push_part_yz_a(psc_particles_c_t *pp)
{
  creal dt = psc.dt;
  creal yl = .5f * dt;
  creal zl = .5f * dt;

  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_c_t *part = psc_particles_c_get_one(pp, n);

    creal root = 1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}

void
genc_push_part_yz_a()
{
  struct psc_genc *genc = psc.c_ctx;
  psc_particles_c_t *pp = &genc->pp;

  static int pr;
  if (!pr) {
    pr = prof_register("genc_part_yz_a", 1., 0, psc.pp.n_part * 12 * sizeof(creal));
  }
  prof_start(pr);
  do_genc_push_part_yz_a(pp);
  prof_stop(pr);
}

