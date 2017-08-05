
#include "psc.h"
#include "psc_particles_cbe.h"

#include <mrc_profile.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

void
particles_cbe_alloc(struct psc_particles *prts, particles_cbe_t *pp, int n_part)
{
  void * m;
  int n_alloced = n_part * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  int ierr = posix_memalign(&m, 16, n_alloced * sizeof(*pp->particles));
  assert(ierr == 0);
  pp->particles = (particle_cbe_t *) m;
}

void particles_cbe_realloc(struct psc_particles *prts, particles_cbe_t *pp, int new_n_part)
{
  if (new_n_part <= psc_mparticles_n_alloced(prts->mprts, prts->p))
    return;
  void * m; 
  int n_alloced = new_n_part * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  free(pp->particles);
  int ierr = posix_memalign(&m, 16, n_alloced * sizeof(*pp->particles));
  assert(ierr == 0);
  pp->particles = (particle_cbe_t *) m;
}

void particles_cbe_free(struct psc_particles *prts, particles_cbe_t *pp)
{
  free(pp->particles);
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, 0);
  pp->particles = NULL;
}

void
psc_mparticles_cbe_get_cbe(struct psc_mparticles *particles, void *_particles_base)
{
  struct psc_mparticles *particles_base = _particles_base;
  *particles = *particles_base;
}

void
psc_mparticles_cbe_put_cbe(struct psc_mparticles *particles, void *particles_base)
{
}

#if 0
static bool __gotten;
#endif

void
psc_mparticles_c_get_cbe(struct psc_mparticles *particles, void *_particles_base)
{
#if 0
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_cbe_get", 1., 0, 0);
  }
  prof_start(pr);


  assert(!__gotten);
  __gotten = true;

  mparticles_c_t *particles_base = _particles_base;
  // With some recent changes, it appears we are allocating/freeing
  // the particles each time. This should encourage us to stop switching
  // between languages/types for the different modules.

  // I cannot currently think of any reason we would offload 
  // the list of patches to the spes, so let's just use calloc 
  // for this one.
  particles->data = calloc(ppsc->nr_patches, sizeof(particles_cbe_t));
  psc_foreach_patch(ppsc, p) {
    particles_c_t *pp_base = psc_mparticles_get_patch_c(particles_base, p);
    particles_cbe_t *pp = psc_mparticles_get_patch_cbe(particles, p);
    pp->n_part = pp_base->n_part;
    // These will be heading to the spes, so we need some memalign lovin'
    void *m;
    int ierr;
    ierr = posix_memalign(&m, 128, pp->n_part * sizeof(*pp->particles));
    pp->particles = (particle_cbe_t *)m;
    assert(ierr == 0);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_c_t *part_base = particles_c_get_one(pp_base,n);
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
    }
  }
#endif
}
      


void
psc_mparticles_c_cbe_put_cbe(struct psc_mparticles *particles, void *_particles_base)
{
#if 0
  assert(__gotten);
  __gotten = false;
  
  mparticles_c_t *particles_base = _particles_base;
  psc_foreach_patch(ppsc, p) {
    particles_c_t *pp_base = psc_mparticles_get_patch_c(particles_base, p);
    particles_cbe_t *pp = psc_mparticles_get_patch_cbe(particles, p);
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++){
      particle_c_t *part_base = particles_c_get_one(pp_base,n);
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
  free(particles->data);
#endif
}


