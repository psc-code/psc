
#include "psc_output_fields_item_private.h"

#include "psc_moments.h"

// ======================================================================

static void
calc_densities(struct psc_output_fields_item *item, mfields_base_t *flds,
	       mparticles_base_t *particles, mfields_c_t *res)
{
  psc_moments_calc_densities(ppsc->moments, flds, particles, res);
}

struct psc_output_fields_item_ops psc_output_fields_item_densities_ops = {
  .name      = "n",
  .nr_comp   = 3,
  .fld_names = { "ne", "ni", "nn" },
  .run       = calc_densities,
};

// ======================================================================

static void
calc_v(struct psc_output_fields_item *item, mfields_base_t *flds,
       mparticles_base_t *particles, mfields_c_t *res)
{
  psc_moments_calc_v(ppsc->moments, flds, particles, res);
}

struct psc_output_fields_item_ops psc_output_fields_item_v_ops = {
  .name      = "v",
  .nr_comp   = 6,
  .fld_names = { "vex", "vey", "vez", "vix", "viy", "viz" },
  .run       = calc_v,
};

// ======================================================================

static void
calc_vv(struct psc_output_fields_item *item, mfields_base_t *flds,
	mparticles_base_t *particles, mfields_c_t *res)
{
  psc_moments_calc_vv(ppsc->moments, flds, particles, res);
}

struct psc_output_fields_item_ops psc_output_fields_item_vv_ops = {
  .name      = "vv",
  .nr_comp   = 6,
  .fld_names = { "vexvex", "veyvey", "vezvez", "vixvix", "viyviy", "vizviz" },
  .run       = calc_vv,
};

