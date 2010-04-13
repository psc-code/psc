
#include "psc.h"

#include <stdlib.h>

static void
genc_particles_from_fortran()
{
  if (!psc.c_part) {
    psc.c_part = calloc(psc.n_part, sizeof(*psc.c_part));
  }
  for (int n = 0; n < psc.n_part; n++) {
    struct f_particle *f_part = &psc.f_part[n];
    struct c_particle *part = &psc.c_part[n];

    part[n].xi  = f_part[n].xi;
    part[n].yi  = f_part[n].yi;
    part[n].zi  = f_part[n].zi;
    part[n].pxi = f_part[n].pxi;
    part[n].pyi = f_part[n].pyi;
    part[n].pzi = f_part[n].pzi;
    part[n].qni = f_part[n].qni;
    part[n].mni = f_part[n].mni;
    part[n].wni = f_part[n].wni;
  }
}

static void
genc_particles_to_fortran()
{
  for (int n = 0; n < psc.n_part; n++) {
    struct f_particle *f_part = &psc.f_part[n];
    struct c_particle *part = &psc.c_part[n];

    f_part[n].xi  = part[n].xi;
    f_part[n].yi  = part[n].yi;
    f_part[n].zi  = part[n].zi;
    f_part[n].pxi = part[n].pxi;
    f_part[n].pyi = part[n].pyi;
    f_part[n].pzi = part[n].pzi;
    f_part[n].qni = part[n].qni;
    f_part[n].mni = part[n].mni;
    f_part[n].wni = part[n].wni;
  }
}

struct psc_ops psc_ops_generic_c = {
  .name = "generic_c",
  .particles_from_fortran = genc_particles_from_fortran,
  .particles_to_fortran   = genc_particles_to_fortran,
  .push_part_yz_a         = genc_push_part_yz_a,
};
