
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
  assert(ppsc->nr_patches == 1);
  static int pr;
  if (!pr) {
    pr = prof_register("fort_out_field", 1., 0, 0);
  }
  prof_start(pr);
  mfields_fortran_t *flds = psc_mfields_fortran_get_from(NE, HZ + 1, flds_base);

  fields_fortran_t *pf = psc_mfields_fortran_get_patch_fortran(flds, 0);
  OUT_field(pf);

  psc_mfields_fortran_put_to(flds, EX, HZ + 1, flds_base);
  prof_stop(pr);
}

// ======================================================================
// psc_output_fields: subclass "fortran"

struct psc_output_fields_ops psc_output_fields_fortran_ops = {
  .name                  = "fortran",
  .run                   = psc_output_fields_fortran_run,
};
