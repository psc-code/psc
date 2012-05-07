
#include "psc_moments_private.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include "psc_moments_1st_cc.c"

// ======================================================================
// psc_moments: subclass "c_1st_cc"

struct psc_moments_ops psc_moments_c_1st_cc_ops = {
  .name                  = "c_1st_cc",
  .calc_densities        = psc_moments_1st_cc_calc_densities,
  .calc_v                = psc_moments_1st_cc_calc_v,
  .calc_vv               = psc_moments_1st_cc_calc_vv,
};

