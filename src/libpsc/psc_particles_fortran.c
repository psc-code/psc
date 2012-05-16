
#include "psc.h"
#include "psc_particles_fortran.h"
#include "psc_particles_c.h"

#include <stdlib.h>

static void
_psc_mparticles_fortran_alloc_patch(mparticles_fortran_t *mp, int p, int n_part)
{
  particles_fortran_t *pp = psc_mparticles_get_patch_fortran(mp, p);
  pp->n_part = n_part;
  pp->n_alloced = n_part * 1.2;
  pp->particles = calloc(pp->n_alloced, sizeof(*pp->particles));
}

static void
_psc_mparticles_fortran_free_patch(mparticles_fortran_t *mp, int p)
{
  particles_fortran_t *pp = psc_mparticles_get_patch_fortran(mp, p);
  free(pp->particles);
  pp->n_alloced = 0;
  pp->particles = NULL;
}

void
particles_fortran_realloc(particles_fortran_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

static void
_psc_mparticles_fortran_copy_to_c(struct psc_mparticles *particles_base,
				  mparticles_c_t *particles, unsigned int flags)
{
  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    pp->n_part = pp_base->n_part;
    assert(pp->n_part <= pp->n_alloced);
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
}

static void
_psc_mparticles_fortran_copy_from_c(struct psc_mparticles *particles_base,
				    mparticles_c_t *particles, unsigned int flags)
{
  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    pp_base->n_part = pp->n_part;
    assert(pp_base->n_part <= pp_base->n_alloced);
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
}

static int
_psc_mparticles_fortran_nr_particles_by_patch(mparticles_fortran_t *mparticles, int p)
{
  return psc_mparticles_get_patch_fortran(mparticles, p)->n_part;
}
									
// ======================================================================
// psc_mparticles: subclass "fortran"
  
static struct mrc_obj_method _psc_mparticles_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",         _psc_mparticles_fortran_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c",       _psc_mparticles_fortran_copy_from_c),
  {}
};

struct psc_mparticles_ops psc_mparticles_fortran_ops = {
  .name                    = "fortran",
  .methods                 = _psc_mparticles_fortran_methods,
  .nr_particles_by_patch   = _psc_mparticles_fortran_nr_particles_by_patch,
  .alloc_patch             = _psc_mparticles_fortran_alloc_patch,
  .free_patch              = _psc_mparticles_fortran_free_patch,
  .size_of_particles_t     = sizeof(particles_fortran_t),
};

