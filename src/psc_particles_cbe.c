
#include "psc.h"
#include "psc_particles_cbe.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#if PARTICLES_BASE == PARTICLES_CBE

static int __arr_size;

void
particles_cbe_alloc(particles_cbe_t *pp, int n_part)
{
  void * m;
  __arr_size = n_part * 1.2;
  //printf("Called to alloc %g particles\n", __arr_size);
 int ierr = posix_memalign(&m, 16, __arr_size * sizeof(*pp->particles));
  assert(ierr == 0);
  pp->particles = (particle_cbe_t *) m;
}

void particles_cbe_realloc(particles_cbe_t *pp, int new_n_part)
{
  if (__arr_size >= new_n_part)
    return;
  void * m; 
  __arr_size = new_n_part * 1.2;
  free(pp->particles);
  int ierr = posix_memalign(&m, 128, __arr_size * sizeof(*pp->particles));
  assert(ierr == 0);
  pp->particles = (particle_cbe_t *) m;
}

void particles_cbe_free(particles_cbe_t *pp)
{
  free(pp->particles);
  pp->particles = NULL;
}

void
particles_cbe_get(mparticles_cbe_t *particles, void *_particles_base)
{
  mparticles_base_t *particles_base = _particles_base;
  *particles = *particles_base;
}

void
particles_cbe_put(mparticles_cbe_t *particles, void *particles_base)
{
}

#else

static particle_cbe_t *__arr;
static int __arr_size;
static bool __gotten;

void
particles_cbe_get(mparticles_cbe_t *particles, void *_particles_base)
{
  assert(!__gotten);
  __gotten = true;

  mparticles_base_t *particles_base = _particles_base;
  // With some recent changes, it appears we are allocating/freeing
  // the particles each time. This should encourage us to stop switching
  // between languages/types for the different modules.

  // I cannot currently think of any reason we would offload 
  // the list of patches to the spes, so let's just use calloc 
  // for this one.
  particles->p = calloc(psc.nr_patches, sizeof(*particles->p));
  foreach_patch(p) {
    particles_base_t *pp_base = &particles_base->p[p];
    particles_cbe_t *pp = &particles->p[p];
    pp->n_part = pp_base->n_part;
    // These will be heading to the spes, so we need some memalign lovin'
    void *m;
    int ierr;
    ierr = posix_memalign(&m, 16, pp->n_part * sizeof(*pp->particles));
    pp->particles = (particle_cbe_t *)m;
    assert(ierr == 0);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_base_t *part_base = particles_base_get_one(pp_base,n);
      particle_cbe_t *part = particles_cbe_get_one(pp,n);

      part->xi  = part_base->xi;
      part->yi  = part_base->yi;
      part->zi  = part_base->zi;
      part->pxi = part_base->pxi;
      part->pyi = part_base->pyi;
      part->pzi = part_base->pzi;
      part->qni = part_base->qni;
      part->mni = part_base->mni;
      part->wni = part_base->wni;
#if PARTICLES_BASE != PARTICLES_C
      part->cni = part_base->cni;
#endif
    }
  }
}
      


void
particles_cbe_put(mparticles_cbe_t *particles, void *_particles_base)
{
  assert(__gotten);
  __gotten = false;
  
  mparticles_base_t *particles_base = _particles_base;
  foreach_patch(p) {
    particles_base_t *pp_base = &particles_base->p[p];
    particles_cbe_t *pp = &particles->p[p];
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++){
      particle_base_t *part_base = particles_base_get_one(pp_base,n);
      particle_cbe_t *part = particles_cbe_get_one(pp,n);

      part_base->xi  = part->xi;
      part_base->yi  = part->yi;
      part_base->zi  = part->zi;
      part_base->pxi = part->pxi;
      part_base->pyi = part->pyi;
      part_base->pzi = part->pzi;
      part_base->qni = part->qni;
      part_base->mni = part->mni;
      part_base->wni = part->wni;
    }
    free(pp->particles);
    pp->particles = NULL;
  }
  free(particles->p);
  particles->p = NULL;
    
}


#endif
