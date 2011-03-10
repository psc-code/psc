
#include "psc_push_fields_private.h"

#include "psc.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_push_fields_fortran_step_a

static void
psc_push_fields_fortran_step_a(struct psc_push_fields *push,
			       mfields_base_t *flds_base)
{
  if (psc.domain.use_pml) {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, MU + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_pml_a", 1., 0, 0);
    }
    prof_start(pr);
    PIC_pml_msa(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, BZ + 1, flds_base);
  } else {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_a", 1., 0, 0);
    }
    prof_start(pr);
    PIC_msa(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, HZ + 1, flds_base);
  }
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_step_b

static void
psc_push_fields_fortran_step_b(struct psc_push_fields *push,
			       mfields_base_t *flds_base)
{
  if (psc.domain.use_pml) {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, MU + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_pml_b", 1., 0, 0);
    }
    prof_start(pr);
    PIC_pml_msb(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, BZ + 1, flds_base);
  } else {
    mfields_fortran_t flds;
    fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
    
    static int pr;
    if (!pr) {
      pr = prof_register("fort_field_b", 1., 0, 0);
    }
    prof_start(pr);
    PIC_msb(&flds.f[0]);
    prof_stop(pr);
    
    fields_fortran_put(&flds, EX, HZ + 1, flds_base);
  }
}

// ======================================================================
// psc_push_fields: subclass "fortran"

struct psc_push_fields_ops psc_push_fields_fortran_ops = {
  .name                  = "fortran",
  .step_a                = psc_push_fields_fortran_step_a,
  .step_b                = psc_push_fields_fortran_step_b,
};
