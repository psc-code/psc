
#include "psc_output_fields_c.h"

#include <mrc_profile.h>

#define to_psc_output_fields_c(out) ((struct psc_output_fields_c *)((out)->obj.subctx))

// ----------------------------------------------------------------------
// psc_output_fields_c_create

static void
psc_output_fields_c_create(struct mrc_obj *obj)
{
  struct psc_output_fields *out = to_psc_output_fields(obj);
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  output_c_create(out_c);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			mfields_base_t *flds_base,
			mparticles_base_t *particles_base)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  output_c_field(out_c, flds_base, particles_base);
}

// ======================================================================
// psc_output_fields: subclass "c"

struct psc_output_fields_ops psc_output_fields_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_output_fields_c),
  .create                = psc_output_fields_c_create,
  .run                   = psc_output_fields_c_run,
};
