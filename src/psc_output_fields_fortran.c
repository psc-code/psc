
#include "psc_output_fields_private.h"

#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_output_fields_fortran_run

static void
psc_output_fields_fortran_run(struct psc_output_fields *out,
			      mfields_base_t *flds_base,
			      mparticles_base_t *particles_base)
{
  assert(psc.nr_patches == 1);
  assert(FIELDS_BASE == FIELDS_FORTRAN);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_out_field", 1., 0, 0);
  }
  prof_start(pr);
  OUT_field();
  prof_stop(pr);
}

// ======================================================================
// psc_output_fields: subclass "fortran"

struct psc_output_fields_ops psc_output_fields_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_output_fields_fortran_run,
};
