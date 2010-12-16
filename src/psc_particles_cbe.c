
#include "psc.h"
#include "psc_particles_cbe.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#if PARTICLES_BASE == PARTICLES_CBE

static int __arr_size;
static int __pad_size;

void
particles_cbe_alloc(particles_cbe_t *pp, int n_part)
{
  void * m;
  __arr_size = n_part * 1.2;
  int ierr = posix_memalign(&m, 128, __arr_size * sizeof(*pp->particles));
  assert(ierr = 0);
  pp->particles = (cbe_particle_t *) m;
}

void particles_cbe_realloc(particles_cbe_t *pp, int new_n_part)
{
  if (__arr_size >= new_n_part)
    return;

  __arr_size = new_n_part * 1.2;
  free(pp->particles);
  int ierr = posix_memalign(&m, 128, __arr_size * sizeof(*pp->particles));
  assert(ierr = 0);
}

void particles_cbe_free(particles_cbe_t *pp)
{
  free(pp->particles);
  pp->particles = NULL;
}

void
particles_cbe_get(particles_cbe_t *pp)
{
  pp->particles = psc.pp.particles;
  pp->n_part = psc.pp.n_part;
}

void
particles_cbe_put(particles_c_t *pp)
{
}

#else

static particle_cbe_t *__arr;
static int __arr_size;
static int __pad_size;
static int __gotten;

void
particles_cbe_get(particles_cbe_t *pp)
{
  void *m;
  int ierr;
  if (psc.pp.n_part > __arr_size) {
    free(__arr);
    __arr = NULL;
  }
  if (!__arr) {
    __arr_size = psc.pp.n_part * 1.2;
    ierr = posix_memalign(&m, 128, __arr_size * sizeof(*__arr));
    __arr = (particle_cbe_t *) m;
  }
  
  assert(!__gotten);
  __gotten = 1;
  pp->particles = __arr;
  pp->null_particles = __pad;
  
  for(int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *f_part = particles_base_get_one(&psc.pp, n);
    particle_cbe_t *part = particles_cbe_get_one(pp,n);

    part->xi  = f_part->xi;
    part->yi  = f_part->yi;
    part->zi  = f_part->zi;
    part->pxi = f_part->pxi;
    part->pyi = f_part->pyi;
    part->pzi = f_part->pzi;
    part->qni = f_part->qni;
    part->mni = f_part->mni;
    part->wni = f_part->wni;
    part->cni = f_part->cni;
  }
}

void
particles_cbe_put(particles_cbe_t *pp)
{
  assert(__gotten);
  __gotten = 0;
  
  for(int n = 0; n < psc.pp.n_part; n++){
    particle_base_t *f_part = particles_base_get_one(&psc.pp, n);
    particle_cbe_t *part = particles_cbe_get_one(pp,n);

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
