
#include "psc.h"
#include "psc_particles_as_fortran.h"
#include "psc_particles_inc.h"
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

// ======================================================================
// conversion to/from "c"

static void
put_particle_c(particle_fortran_t *prt, int n, struct psc_particles *prts_c)
{
  particle_c_t *prt_c = particles_c_get_one(prts_c, n);
  
  prt_c->xi      = prt->xi;
  prt_c->yi      = prt->yi;
  prt_c->zi      = prt->zi;
  prt_c->pxi     = prt->pxi;
  prt_c->pyi     = prt->pyi;
  prt_c->pzi     = prt->pzi;
  prt_c->qni     = prt->qni;
  prt_c->wni     = prt->wni;
  prt_c->mni     = prt->mni;
}

static void
get_particle_c(particle_fortran_t *prt, int n, struct psc_particles *prts_c)
{
  particle_c_t *prt_c = particles_c_get_one(prts_c, n);

  prt->xi      = prt_c->xi;
  prt->yi      = prt_c->yi;
  prt->zi      = prt_c->zi;
  prt->pxi     = prt_c->pxi;
  prt->pyi     = prt_c->pyi;
  prt->pzi     = prt_c->pzi;
  prt->qni     = prt_c->qni;
  prt->wni     = prt_c->wni;
  prt->mni     = prt_c->mni;
}

static void
psc_mparticles_fortran_copy_to_c(int p, struct psc_mparticles *mprts,
				 struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts, mprts_c, flags, put_particle_c);
}

static void
psc_mparticles_fortran_copy_from_c(int p, struct psc_mparticles *mprts,
				   struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts, mprts_c, flags, get_particle_c);
}

// ======================================================================
// conversion to/from "double"

static inline void
calc_vxi(particle_fortran_real_t vxi[3], particle_fortran_t *part)
{
  particle_fortran_real_t root =
    1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static void
put_particle_double(particle_fortran_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_fortran_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_double_t *prt_dbl = particles_double_get_one(prts_dbl, n);

  particle_double_real_t vxi[3];
  calc_vxi(vxi, prt);
  
  prt_dbl->xi      = prt->xi + dth[0] * vxi[0];
  prt_dbl->yi      = prt->yi + dth[1] * vxi[1];
  prt_dbl->zi      = prt->zi + dth[2] * vxi[2];
  prt_dbl->pxi     = prt->pxi;
  prt_dbl->pyi     = prt->pyi;
  prt_dbl->pzi     = prt->pzi;
  prt_dbl->qni_wni = prt->qni * prt->wni;;
  prt_dbl->kind    = prt->qni > 0 ? 1 : 0;
}

static void
get_particle_double(particle_fortran_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_fortran_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_double_t *prt_dbl = particles_double_get_one(prts_dbl, n);

  particle_fortran_real_t qni = ppsc->kinds[prt_dbl->kind].q;
  particle_fortran_real_t mni = ppsc->kinds[prt_dbl->kind].m;
  particle_fortran_real_t wni = prt_dbl->qni_wni / qni;
  
  prt->xi  = prt_dbl->xi;
  prt->yi  = prt_dbl->yi;
  prt->zi  = prt_dbl->zi;
  prt->pxi = prt_dbl->pxi;
  prt->pyi = prt_dbl->pyi;
  prt->pzi = prt_dbl->pzi;
  prt->qni = qni;
  prt->mni = mni;
  prt->wni = wni;

  particle_fortran_real_t vxi[3];
  calc_vxi(vxi, prt);
  prt->xi -= dth[0] * vxi[0];
  prt->yi -= dth[1] * vxi[1];
  prt->zi -= dth[2] * vxi[2];
}

static void
psc_mparticles_fortran_copy_to_double(int p, struct psc_mparticles *mprts_fortran,
				      struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts_fortran, mprts_dbl, flags, put_particle_double);
}

static void
psc_mparticles_fortran_copy_from_double(int p, struct psc_mparticles *mprts_fortran,
					struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts_fortran, mprts_dbl, flags, get_particle_double);
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

