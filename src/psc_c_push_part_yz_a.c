
#include "psc.h"

#include <math.h>
#include <stdlib.h>

// test w/ C data structure, single precision 

void
psc_c_push_part_yz_a()
{
  real dt = psc.dt;
  real yl = .5 * dt;
  real zl = .5 * dt;

  for (int n = 0; n < psc.n_part; n++) {
    struct c_particle *part = &psc.c_part[n];

    real root = 1. / sqrt(1. + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    real vyi = part->pyi * root;
    real vzi = part->pzi * root;

    part->yi += vyi * yl;
    part->zi += vzi * zl;
  }
}

void
psc_push_part_yz_a_c()
{
  psc_particles_from_fortran();
  psc_c_push_part_yz_a();
  psc_particles_to_fortran();
}
