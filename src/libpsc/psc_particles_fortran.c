
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
_psc_mparticles_fortran_get_fortran(void *_particles_base)
{
  return _particles_base;
}

void
_psc_mparticles_fortran_put_fortran(mparticles_fortran_t *particles, void *_particles_base)
{
}

void
particles_fortran_realloc(particles_fortran_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

static bool __gotten;

static mparticles_c_t *
_psc_mparticles_fortran_get_c(void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;
    
  mparticles_fortran_t *particles_base = _particles_base;

  int *nr_particles_by_patch = malloc(particles_base->nr_patches * sizeof(int));
  for (int p = 0; p < particles_base->nr_patches; p++) {
    nr_particles_by_patch[p] = psc_mparticles_fortran_nr_particles_by_patch(particles_base, p);
  }
  struct mrc_domain *domain = particles_base->domain;
  mparticles_c_t *particles = psc_mparticles_c_create(mrc_domain_comm(domain));
  psc_mparticles_c_set_domain_nr_particles(particles, domain, nr_particles_by_patch);
  psc_mparticles_c_setup(particles);
  free(nr_particles_by_patch);

  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    pp->n_part = pp_base->n_part;
    pp->particles = calloc(pp->n_part, sizeof(particles_fortran_t));
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *part_base = particles_fortran_get_one(pp_base, n);
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

  prof_stop(pr);
  return particles;
}

static void
_psc_mparticles_fortran_put_c(mparticles_c_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_fortran_t *particles_base = _particles_base;
  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *part_base = particles_fortran_get_one(pp_base, n);
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
  }
  psc_mparticles_c_destroy(particles);

  prof_stop(pr);
}

static void
_psc_mparticles_fortran_set_domain_nr_particles(mparticles_fortran_t *mparticles,
						struct mrc_domain *domain,
						int *nr_particles_by_patch)
{
  mparticles->domain = domain;
  mrc_domain_get_patches(domain, &mparticles->nr_patches);

  mparticles->data = calloc(mparticles->nr_patches,
			    sizeof(*mparticles->data));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_fortran_alloc(psc_mparticles_get_patch_fortran(mparticles, p),
			    nr_particles_by_patch[p]);
  }
}
									
static void
_psc_mparticles_fortran_destroy(mparticles_fortran_t *mparticles)
{
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_fortran_free(psc_mparticles_get_patch_fortran(mparticles, p));
  }
  free(mparticles->data);
}

static int
_psc_mparticles_fortran_nr_particles_by_patch(mparticles_fortran_t *mparticles, int p)
{
  return psc_mparticles_get_patch_fortran(mparticles, p)->n_part;
}
									
// ======================================================================
// psc_mparticles: subclass "fortran"
  
struct psc_mparticles_fortran_ops psc_mparticles_fortran_ops = {
  .name                    = "fortran",
  .set_domain_nr_particles = _psc_mparticles_fortran_set_domain_nr_particles,
  .nr_particles_by_patch   = _psc_mparticles_fortran_nr_particles_by_patch,
  .get_c                   = _psc_mparticles_fortran_get_c,
  .put_c                   = _psc_mparticles_fortran_put_c,
  .get_fortran             = _psc_mparticles_fortran_get_fortran,
  .put_fortran             = _psc_mparticles_fortran_put_fortran,
};

static void
psc_mparticles_fortran_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles_fortran, &psc_mparticles_fortran_ops);
}

struct mrc_class_psc_mparticles_fortran mrc_class_psc_mparticles_fortran = {
  .name             = "psc_mparticles_fortran",
  .size             = sizeof(struct psc_mparticles),
  .init             = psc_mparticles_fortran_init,
  .destroy          = _psc_mparticles_fortran_destroy,
};

