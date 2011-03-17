
#include "psc.h"
#include "psc_particles_c.h"

#include <stdlib.h>
#include <assert.h>

#if PARTICLES_BASE == PARTICLES_C

void
particles_c_alloc(particles_c_t *pp, int n_part)
{
  pp->n_alloced = n_part * 1.2;
  pp->particles = calloc(pp->n_alloced, sizeof(*pp->particles));
}

void
particles_c_realloc(particles_c_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

void
particles_c_free(particles_c_t *pp)
{
  free(pp->particles);
  pp->n_alloced = 0;
  pp->particles = NULL;
}

void
particles_c_get(mparticles_c_t *particles, void *_particles_base)
{
  mparticles_base_t *particles_base = _particles_base;
  *particles = *particles_base;
}

void
particles_c_put(mparticles_c_t *particles, void *particles_base)
{
}

#else

static bool __gotten;

void
particles_c_get(mparticles_c_t *particles, void *_particles_base)
{
  assert(!__gotten);
  __gotten = true;
    
  mparticles_base_t *particles_base = _particles_base;

  particles->p = calloc(psc.nr_patches, sizeof(*particles->p));
  foreach_patch(p) {
    particles_base_t *pp_base = &particles_base->p[p];
    particles_c_t *pp = &particles->p[p];
    pp->n_part = pp_base->n_part;
    pp->particles = calloc(pp->n_part, sizeof(*pp->particles));
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_base_t *part_base = particles_base_get_one(pp_base, n);
      particle_c_t *part = particles_c_get_one(pp, n);
      
      part->xi  = part_base->xi;
      part->yi  = part_base->yi;
      part->zi  = part_base->zi;
      part->pxi = part_base->pxi;
      part->pyi = part_base->pyi;
      part->pzi = part_base->pzi;
      part->qni = part_base->qni;
      part->mni = part_base->mni;
      part->wni = part_base->wni;
    }
  }
}

void
particles_c_put(mparticles_c_t *particles, void *_particles_base)
{
  assert(__gotten);
  __gotten = false;

  mparticles_base_t *particles_base = _particles_base;
  foreach_patch(p) {
    particles_base_t *pp_base = &particles_base->p[p];
    particles_c_t *pp = &particles->p[p];
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_base_t *part_base = particles_base_get_one(pp_base, n);
      particle_c_t *part = particles_c_get_one(pp, n);
      
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
  }
  free(particles->p);
  particles->p = NULL;
}

#endif

