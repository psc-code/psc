
#include "psc.h"
#include "psc_particles_fortran.h"
#include "psc_particles_c.h"
#include "psc_particles_double.h"

#include <stdlib.h>

// ======================================================================
// psc_particles "fortran"

static void
psc_particles_fortran_setup(struct psc_particles *prts)
{
  struct psc_particles_fortran *fort = psc_particles_fortran(prts);

  prts->n_alloced = psc_particles_size(prts) * 1.2;
  fort->particles = calloc(prts->n_alloced, sizeof(*fort->particles));
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

  if (new_n_part <= prts->n_alloced)
    return;

  prts->n_alloced = new_n_part * 1.2;
  fort->particles = realloc(fort->particles, prts->n_alloced * sizeof(*fort->particles));
}

static void
psc_mparticles_fortran_copy_to_c(int p, struct psc_mparticles *mprts_base,
				 struct psc_mparticles *mprts_c, unsigned int flags)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(mprts_base, p);
  struct psc_particles *prts_c = psc_mparticles_get_patch(mprts_c, p);
  int n_prts = psc_particles_size(prts_base);
  psc_particles_resize(prts_c, n_prts);
  for (int n = 0; n < n_prts; n++) {
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
psc_mparticles_fortran_copy_from_c(int p, struct psc_mparticles *mprts_base,
				   struct psc_mparticles *mprts_c, unsigned int flags)
{
  struct psc_particles *prts_base = psc_mparticles_get_patch(mprts_base, p);
  struct psc_particles *prts_c = psc_mparticles_get_patch(mprts_c, p);
  int n_prts = psc_particles_size(prts_c);
  psc_particles_resize(prts_base, n_prts);
  for (int n = 0; n < n_prts; n++) {
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

static inline void
calc_vxi(particle_double_real_t vxi[3], particle_double_t *part)
{
  particle_double_real_t root =
    1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static void
psc_mparticles_fortran_copy_to_double(int p, struct psc_mparticles *mprts_fortran,
				      struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  struct psc_particles *prts_fortran = psc_mparticles_get_patch(mprts_fortran, p);
  struct psc_particles *prts_dbl = psc_mparticles_get_patch(mprts_dbl, p);
  particle_double_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }

  int n_prts = psc_particles_size(prts_fortran);
  psc_particles_resize(prts_dbl, n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_double_t *part_dbl = particles_double_get_one(prts_dbl, n);
    particle_fortran_t *part_fortran = particles_fortran_get_one(prts_fortran, n);
    
    particle_double_real_t qni_wni;
    if (part_fortran->qni != 0.) {
      qni_wni = part_fortran->qni * part_fortran->wni;
    } else {
      qni_wni = part_fortran->wni;
    }
    
    part_dbl->xi          = part_fortran->xi;
    part_dbl->yi          = part_fortran->yi;
    part_dbl->zi          = part_fortran->zi;
    part_dbl->pxi         = part_fortran->pxi;
    part_dbl->pyi         = part_fortran->pyi;
    part_dbl->pzi         = part_fortran->pzi;
    part_dbl->qni_wni     = qni_wni;
    part_dbl->kind        = (qni_wni > 0) ? 1 : 0;

    particle_double_real_t vxi[3];
    calc_vxi(vxi, part_dbl);
    part_dbl->xi += dth[0] * vxi[0];
    part_dbl->yi += dth[1] * vxi[1];
    part_dbl->zi += dth[2] * vxi[2];
  }
}

static void
psc_mparticles_fortran_copy_from_double(int p, struct psc_mparticles *mprts_fortran,
					struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  struct psc_particles *prts_fortran = psc_mparticles_get_patch(mprts_fortran, p);
  struct psc_particles *prts_dbl = psc_mparticles_get_patch(mprts_dbl, p);
  particle_double_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }

  int n_prts = psc_particles_size(prts_dbl);
  psc_particles_resize(prts_fortran, n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_double_t *part_dbl = particles_double_get_one(prts_dbl, n);
    particle_fortran_t *part_fortran = particles_fortran_get_one(prts_fortran, n);
    
    particle_fortran_real_t qni = ppsc->kinds[part_dbl->kind].q;
    particle_fortran_real_t mni = ppsc->kinds[part_dbl->kind].m;
    particle_fortran_real_t wni = part_dbl->qni_wni / qni;
    
    particle_double_real_t vxi[3];
    calc_vxi(vxi, part_dbl);
    part_fortran->xi  = part_dbl->xi - dth[0] * vxi[0];
    part_fortran->yi  = part_dbl->yi - dth[1] * vxi[1];
    part_fortran->zi  = part_dbl->zi - dth[2] * vxi[2];
    part_fortran->pxi = part_dbl->pxi;
    part_fortran->pyi = part_dbl->pyi;
    part_fortran->pzi = part_dbl->pzi;
    part_fortran->qni = qni;
    part_fortran->mni = mni;
    part_fortran->wni = wni;
  }
}

// ======================================================================
// psc_particles: subclass "fortran"

struct psc_particles_ops psc_particles_fortran_ops = {
  .name                    = "fortran",
  .size                    = sizeof(struct psc_particles_fortran),
  .setup                   = psc_particles_fortran_setup,
  .destroy                 = psc_particles_fortran_destroy,
};

// ======================================================================
// psc_mparticles: subclass "fortran"
  
static struct mrc_obj_method psc_particles_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mparticles_fortran_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mparticles_fortran_copy_from_c),
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_fortran_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_fortran_copy_from_double),
  {}
};

struct psc_mparticles_ops psc_mparticles_fortran_ops = {
  .name                    = "fortran",
  .methods                 = psc_particles_fortran_methods,
};

