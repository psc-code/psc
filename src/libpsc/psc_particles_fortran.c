
#include "psc.h"
#include "psc_particles_fortran.h"
#include "psc_particles_c.h"

#include <stdlib.h>

static void *
_psc_mparticles_fortran_alloc_patch(int p, int n_part, unsigned int flags)
{
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME!
  struct psc_particles *prts = psc_particles_create(comm);
  psc_particles_set_type(prts, "fortran");
  prts->n_part = n_part;
  psc_particles_setup(prts);
  return prts;
}

static void
_psc_mparticles_fortran_free_patch(int p, void *_pp)
{
  struct psc_particles *prts = _pp;
  psc_particles_destroy(prts);
}

// ======================================================================
// psc_particles "fortran"

static void
psc_particles_fortran_setup(struct psc_particles *prts)
{
  struct psc_particles_fortran *fort = psc_particles_fortran(prts);

  fort->n_alloced = prts->n_part * 1.2;
  fort->particles = calloc(fort->n_alloced, sizeof(*fort->particles));
}

static void
psc_particles_fortran_destroy(struct psc_particles *prts)
{
  struct psc_particles_fortran *fort = psc_particles_fortran(prts);

  free(fort->particles);
}

void
particles_fortran_realloc(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_fortran *fort = psc_particles_fortran(prts);

  if (new_n_part <= fort->n_alloced)
    return;

  fort->n_alloced = new_n_part * 1.2;
  fort->particles = realloc(fort->particles, fort->n_alloced * sizeof(*fort->particles));
}

static void
_psc_mparticles_fortran_copy_to_c(int p, struct psc_mparticles *particles_base,
				  mparticles_c_t *particles, unsigned int flags)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, p);
  struct psc_particles *prts_c = psc_mparticles_get_patch(particles, p);
  struct psc_particles_c *c = psc_particles_c(prts_c);
  prts_c->n_part = prts_base->n_part;
  assert(prts_c->n_part <= c->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_fortran_t *part_base = particles_fortran_get_one(prts_base, n);
    particle_c_t *part = particles_c_get_one(prts_c, n);
    
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

static void
_psc_mparticles_fortran_copy_from_c(int p, struct psc_mparticles *particles_base,
				    mparticles_c_t *particles, unsigned int flags)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_base, p);
  struct psc_particles_fortran *fort = psc_particles_fortran(prts_base);
  struct psc_particles *prts_c = psc_mparticles_get_patch(particles, p);
  prts_base->n_part = prts_c->n_part;
  assert(prts_base->n_part <= fort->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_fortran_t *part_base = particles_fortran_get_one(prts_base, n);
    particle_c_t *part = particles_c_get_one(prts_c, n);
    
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

static int
_psc_mparticles_fortran_nr_particles_by_patch(mparticles_fortran_t *mparticles, int p)
{
  return psc_mparticles_get_patch(mparticles, p)->n_part;
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
};

// ======================================================================
// psc_particles: subclass "fortran"

struct psc_particles_ops psc_particles_fortran_ops = {
  .name                    = "fortran",
  .size                    = sizeof(struct psc_particles_fortran),
  .setup                   = psc_particles_fortran_setup,
  .destroy                 = psc_particles_fortran_destroy,
};
