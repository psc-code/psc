
#include "psc_moments_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_moments_fortran_calc_densities

static void
psc_moments_fortran_calc_densities(struct psc_moments *out,
				   mfields_base_t *flds_base,
				   mparticles_base_t *particles_base,
				   mfields_c_t *res)
{
  assert(0);
#if 0
  assert(ppsc->nr_patches == 1);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_densities", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  psc_mparticles_fortran_get_from(&particles, &particles_base);
  mfields_fortran_t flds_fortran;
  psc_mfields_fortran_get_from(&flds_fortran, 0, 0, res);

  CALC_densities(&particles.p[0], &flds_fortran.f[0]);

  psc_mparticles_fortran_put_to(&particles, &ppsc->particles);
  psc_mfields_fortran_put_to(&flds_fortran, NE, NE + 3, res);

  prof_stop(pr);
#endif
}

// ======================================================================
// psc_moments: subclass "fortran"

struct psc_moments_ops psc_moments_fortran_ops = {
  .name                  = "fortran",
  .calc_densities        = psc_moments_fortran_calc_densities,
};
