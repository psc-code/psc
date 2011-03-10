
#include "psc_output_fields_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_output_fields_c_create

static void
psc_output_fields_c_create(struct mrc_obj *obj)
{
  struct psc_output_fields *out = to_psc_output_fields(obj);
  output_c_create();
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			mfields_base_t *flds_base,
			mparticles_base_t *particles_base)
{
  output_c_field(flds_base, particles_base);
}

// ======================================================================
// psc_output_fields: subclass "c"

struct psc_output_fields_ops psc_output_fields_c_ops = {
  .name                  = "c",
  .create                = psc_output_fields_c_create,
  .run                   = psc_output_fields_c_run,
};
