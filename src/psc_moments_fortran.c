
#include "psc_moments_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_moments_fortran_calc_densities

static void
psc_moments_fortran_calc_densities(struct psc_moments *out,
				   mfields_base_t *flds_base,
				   mparticles_base_t *particles_base,
				   mfields_base_t *res)
{
  assert(psc.nr_patches == 1);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_densities", 1., 0, 0);
  }
  prof_start(pr);

  mparticles_fortran_t particles;
  particles_fortran_get(&particles, &particles_base);
  mfields_fortran_t flds_fortran;
  fields_fortran_get_from(&flds_fortran, 0, 0, res, 0);

  CALC_densities(&particles.p[0], &flds_fortran.f[0]);

  particles_fortran_put(&particles, &psc.particles);
  fields_fortran_put_to(&flds_fortran, NE, NE + 3, res, 0);

  prof_stop(pr);
}

// ======================================================================
// psc_moments: subclass "fortran"

struct psc_moments_ops psc_moments_fortran_ops = {
  .name                  = "fortran",
  .calc_densities        = psc_moments_fortran_calc_densities,
};
