
#include "psc.h"
#include "psc_particles_cuda2.h"

// for conversions
#include "psc_particles_single.h"

#include <stdlib.h>

// ======================================================================
// psc_particles "cuda2"

// ----------------------------------------------------------------------
// psc_particles_cuda2_setup

static void
psc_particles_cuda2_setup(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  sub->n_alloced = prts->n_part * 1.2;
  sub->particles = calloc(sub->n_alloced, sizeof(*sub->particles));
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_destroy

static void
psc_particles_cuda2_destroy(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  free(sub->particles);
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_to_single

static void
psc_particles_cuda2_copy_to_single(struct psc_particles *prts_base,
				    struct psc_particles *prts, unsigned int flags)
{
  prts->n_part = prts_base->n_part;
  assert(prts->n_part <= psc_particles_single(prts)->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_cuda2_t *part_base = particles_cuda2_get_one(prts_base, n);
    particle_single_t *part = particles_single_get_one(prts, n);
    
    part->xi      = part_base->xi;
    part->yi      = part_base->yi;
    part->zi      = part_base->zi;
    part->pxi     = part_base->pxi;
    part->pyi     = part_base->pyi;
    part->pzi     = part_base->pzi;
    part->qni_wni = part_base->qni_wni;
    part->kind    = part_base->kind;
  }
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_from_single

static void
psc_particles_cuda2_copy_from_single(struct psc_particles *prts_base,
				      struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts_base);
  prts_base->n_part = prts->n_part;
  assert(prts_base->n_part <= sub->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_cuda2_t *part_base = particles_cuda2_get_one(prts_base, n);
    particle_single_t *part = particles_single_get_one(prts, n);

    part_base->xi      = part->xi;
    part_base->yi      = part->yi;
    part_base->zi      = part->zi;
    part_base->pxi     = part->pxi;
    part_base->pyi     = part->pyi;
    part_base->pzi     = part->pzi;
    part_base->qni_wni = part->qni_wni;
    part_base->kind    = part->kind;
  }
}

// ======================================================================
// psc_particles: subclass "cuda2"

static struct mrc_obj_method psc_particles_cuda2_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_cuda2_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_cuda2_copy_from_single),
  {}
};

struct psc_particles_ops psc_particles_cuda2_ops = {
  .name                    = "cuda2",
  .size                    = sizeof(struct psc_particles_cuda2),
  .methods                 = psc_particles_cuda2_methods,
  .setup                   = psc_particles_cuda2_setup,
  .destroy                 = psc_particles_cuda2_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_cuda2_read,
  .write                   = psc_particles_cuda2_write,
#endif
  .reorder                 = psc_particles_cuda2_reorder,
#endif
};

// ======================================================================
// psc_mparticles: subclass "cuda2"
  
struct psc_mparticles_ops psc_mparticles_cuda2_ops = {
  .name                    = "cuda2",
};

