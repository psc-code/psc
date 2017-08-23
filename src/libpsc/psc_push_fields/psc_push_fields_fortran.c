
#include "psc_push_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_E

static void
psc_push_fields_fortran_push_E(struct psc_push_fields *push,
			       struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", JXI, HZ + 1);
  PIC_msa_e(flds);
  psc_fields_put_as(flds, flds_base, EX, HZ + 1);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_H

static void
psc_push_fields_fortran_push_H(struct psc_push_fields *push,
			       struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", JXI, HZ + 1);
  PIC_msa_h(flds);
  psc_fields_put_as(flds, flds_base, EX, HZ + 1);
}

// ======================================================================
// psc_push_fields: subclass "fortran"

struct psc_push_fields_ops psc_push_fields_fortran_ops = {
  .name                  = "fortran",
  .push_E                = psc_push_fields_fortran_push_E,
  .push_H                = psc_push_fields_fortran_push_H,
};
