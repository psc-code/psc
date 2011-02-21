
#include "psc.h"
#include "psc_particles_c.h"

#include <stdlib.h>
#include <assert.h>

#if PARTICLES_BASE == PARTICLES_C

static int __arr_size;

void
particles_c_alloc(particles_c_t *pp, int n_part)
{
  __arr_size = n_part * 1.2;
  pp->particles = calloc(__arr_size, sizeof(*pp->particles));
}

void
particles_c_realloc(particles_c_t *pp, int new_n_part)
{
  if (__arr_size <= new_n_part)
    return;

  __arr_size = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, __arr_size * sizeof(*pp->particles));
}

void
particles_c_free(particles_c_t *pp)
{
  free(pp->particles);
  pp->particles = NULL;
}

void
particles_c_get(particles_c_t *pp, struct psc_mparticles *particles_base)
{
  *pp = particles_base->p[0];
}

void
particles_c_put(particles_c_t *pp, struct psc_mparticles *particles_base)
{
}

#else

static particle_c_t *__arr;
static int __arr_size;
static int __gotten;

void
particles_c_get(particles_c_t *pp, struct psc_mparticles *particles_base)
{
  particles_base_t *pp_base = &particles_base->p[0];
  if (pp_base->n_part > __arr_size) {
    free(__arr);
    __arr = NULL;
  }
  if (!__arr) {
    __arr_size = pp_base->n_part * 1.2;
    __arr = calloc(__arr_size, sizeof(*__arr));
  }
  assert(!__gotten);
  __gotten = 1;

  pp->particles = __arr;
  for (int n = 0; n < pp_base->n_part; n++) {
    particle_base_t *f_part = particles_base_get_one(pp_base, n);
    particle_c_t *part = particles_c_get_one(pp, n);

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
particles_c_put(particles_c_t *pp, struct psc_mparticles *particles_base)
{
  assert(__gotten);
  __gotten = 0;

  particles_base_t *pp_base = &particles_base->p[0];
  for (int n = 0; n < pp_base->n_part; n++) {
    particle_base_t *f_part = particles_base_get_one(pp_base, n);
    particle_c_t *part = particles_c_get_one(pp, n);

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

  pp->particles = NULL;
}

#endif

