
#include "psc.h"

#include <math.h>

void
genf_push_part_yz_a()
{
  f_real dt = psc.dt;
  f_real yl = .5 * dt;
  f_real zl = .5 * dt;

  for (int n = 0; n < psc.n_part; n++) {
    struct f_particle *part = &psc.f_part[n];

    f_real root = 1. / sqrt(1. + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    f_real vyi = part->pyi * root;
    f_real vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}
