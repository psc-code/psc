
#include "psc_moments_private.h"

#include "psc_glue.h"

// ----------------------------------------------------------------------
// psc_moments_fortran_calc_densities

static void
psc_moments_fortran_calc_densities(struct psc_moments *out,
				   mfields_base_t *flds_base,
				   mparticles_base_t *particles_base,
				   mfields_c_t *res)
{
  assert(ppsc->nr_patches == 1);

  mparticles_fortran_t *particles = psc_mparticles_get_fortran(particles_base, 0);
  mfields_fortran_t *flds_fortran = psc_mfields_get_fortran(res, 0, 0);

  CALC_densities(psc_mparticles_get_patch_fortran(particles, 0),
		 psc_mfields_get_patch_fortran(flds_fortran, 0));

  psc_mparticles_put_fortran(particles, particles_base, MP_DONT_COPY);
  psc_mfields_put_fortran(flds_fortran, res, 0, 3);
}

// ======================================================================
// psc_moments: subclass "fortran"

struct psc_moments_ops psc_moments_fortran_ops = {
  .name                  = "fortran",
  .calc_densities        = psc_moments_fortran_calc_densities,
};
