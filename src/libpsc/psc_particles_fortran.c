
#include "psc.h"
#include "psc_particles_fortran.h"

#include <mrc_profile.h>
#include <stdlib.h>

void
particles_fortran_alloc(particles_fortran_t *pp, int n_part)
{
  pp->n_alloced = n_part * 1.2;
  pp->particles = calloc(pp->n_alloced, sizeof(*pp->particles));
}

void
particles_fortran_free(particles_fortran_t *pp)
{
  free(pp->particles);
  pp->n_alloced = 0;
  pp->particles = NULL;
}

mparticles_fortran_t *
psc_mparticles_fortran_get_fortran(void *_particles_base)
{
  return _particles_base;
}

void
psc_mparticles_fortran_put_fortran(mparticles_fortran_t *particles, void *_particles_base)
{
}

static bool __gotten;

mparticles_fortran_t *
psc_mparticles_c_get_fortran(void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("part_fortran_get", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;

  mparticles_c_t *particles_base = _particles_base;

  mparticles_fortran_t *particles = calloc(1, sizeof(*particles));
  particles->data = calloc(ppsc->nr_patches, sizeof(*particles->data));
  psc_foreach_patch(ppsc, p) {
    particles_c_t *pp_base = psc_mparticles_get_patch_c(particles_base, p);
    particles_fortran_t *pp = psc_mparticles_get_patch_fortran(particles, p);
    pp->n_part = pp_base->n_part;
    pp->particles = calloc(pp->n_part, sizeof(*pp->particles));

    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *f_part = particles_fortran_get_one(pp, n);
      particle_c_t *part = particles_c_get_one(pp_base, n);

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

  prof_stop(pr);
  return particles;
}

void
psc_mparticles_c_put_fortran(mparticles_fortran_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("part_fortran_put", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_c_t *particles_base = _particles_base;
  psc_foreach_patch(ppsc, p) {
    particles_c_t *pp_base = psc_mparticles_get_patch_c(particles_base, p);
    particles_fortran_t *pp = psc_mparticles_get_patch_fortran(particles, p);
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *f_part = &pp->particles[n];
      particle_c_t *part = particles_c_get_one(pp_base, n);
      
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

    free(pp->particles);
  }
  free(particles->data);
  free(particles);

  prof_stop(pr);
}

void
particles_fortran_realloc(particles_fortran_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

