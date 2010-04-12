
#include "psc.h"

#include <math.h>

void
psc_push_part_yz_a()
{
  real dt = psc.dt;
  real yl = .5 * dt;
  real zl = .5 * dt;

  for (int n = 0; n < psc.n_part; n++) {
    struct f_particle *part = &psc.f_part[n];

    real root = 1. / sqrt(1. + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    real vyi = part->pyi * root;
    real vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}
