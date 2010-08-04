
#include "psc.h"
#include "psc_particles_c.h"

#include <stdlib.h>
#include <assert.h>

static particle_c_t *__arr;
static int __arr_size;

void
psc_particles_c_get(psc_particles_c_t *pp)
{
  if (psc.pp.n_part > __arr_size) {
    free(__arr);
    __arr = NULL;
  }
  if (!__arr) {
    __arr_size = psc.pp.n_part * 1.2;
    __arr = calloc(__arr_size, sizeof(*__arr));
  }

  pp->particles = __arr;
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
    particle_c_t *part = psc_particles_c_get_one(pp, n);

    part->xi  = f_part->xi;
    part->yi  = f_part->yi;
    part->zi  = f_part->zi;
    part->pxi = f_part->pxi;
    part->pyi = f_part->pyi;
    part->pzi = f_part->pzi;
    part->qni = f_part->qni;
    part->mni = f_part->mni;
    part->wni = f_part->wni;
  }
}

void
psc_particles_c_put(psc_particles_c_t *pp)
{
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
    particle_c_t *part = psc_particles_c_get_one(pp, n);

    f_part->xi  = part->xi;
    f_part->yi  = part->yi;
    f_part->zi  = part->zi;
    f_part->pxi = part->pxi;
    f_part->pyi = part->pyi;
    f_part->pzi = part->pzi;
    f_part->qni = part->qni;
    f_part->mni = part->mni;
    f_part->wni = part->wni;
  }
}

