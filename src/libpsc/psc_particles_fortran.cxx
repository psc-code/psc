
#include "psc.h"
#include "psc_particles_as_fortran.h"
#include "psc_particles_inc.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_mparticles: subclass "fortran"
  
// ----------------------------------------------------------------------
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
convert_to_double(particle_fortran_t *prt, int n, struct psc_mparticles *mprts_dbl, int p)
{
  particle_fortran_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_double_t *prt_dbl = &mparticles_double_t(mprts_dbl)[p][n];

  particle_double_real_t vxi[3];
  calc_vxi(vxi, prt);
  
  prt_dbl->xi      = prt->xi + dth[0] * vxi[0];
  prt_dbl->yi      = prt->yi + dth[1] * vxi[1];
  prt_dbl->zi      = prt->zi + dth[2] * vxi[2];
  prt_dbl->pxi     = prt->pxi;
  prt_dbl->pyi     = prt->pyi;
  prt_dbl->pzi     = prt->pzi;
  prt_dbl->qni_wni = prt->qni * prt->wni;;
  prt_dbl->kind_   = prt->qni > 0 ? 1 : 0;
}

static void
convert_from_double(particle_fortran_t *prt, int n, struct psc_mparticles *mprts_dbl, int p)
{
  particle_fortran_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_double_t *prt_dbl = &mparticles_double_t(mprts_dbl)[p][n];

  particle_fortran_real_t qni = ppsc->kinds[prt_dbl->kind_].q;
  particle_fortran_real_t mni = ppsc->kinds[prt_dbl->kind_].m;
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
psc_mparticles_fortran_copy_to_double(struct psc_mparticles *mprts_fortran,
				      struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_to(mprts_fortran, mprts_dbl, flags, convert_to_double);
}

static void
psc_mparticles_fortran_copy_from_double(struct psc_mparticles *mprts_fortran,
					struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_from(mprts_fortran, mprts_dbl, flags, convert_from_double);
}

static struct mrc_obj_method psc_mparticles_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_fortran_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_fortran_copy_from_double),
  {}
};

#include "psc_particles_common.cxx"

