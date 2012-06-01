
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

// ----------------------------------------------------------------------
// psc_push_fields_fortran_pml_a

static void
psc_push_fields_fortran_pml_a(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  assert(ppsc->nr_patches == 1);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", JXI, MU + 1);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_field_pml_a", 1., 0, 0);
  }
  prof_start(pr);
  PIC_pml_msa(flds);
  prof_stop(pr);
  
  psc_fields_put_as(flds, flds_base, EX, BZ + 1);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_pml_b

static void
psc_push_fields_fortran_pml_b(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  assert(ppsc->nr_patches == 1);
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", JXI, MU + 1);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_field_pml_b", 1., 0, 0);
  }
  prof_start(pr);
  PIC_pml_msb(flds);
  prof_stop(pr);
  
  psc_fields_put_as(flds, flds_base, EX, BZ + 1);
}

// ======================================================================
// psc_push_fields: subclass "fortran"

struct psc_push_fields_ops psc_push_fields_fortran_ops = {
  .name                  = "fortran",
  .push_E                = psc_push_fields_fortran_push_E,
  .push_H                = psc_push_fields_fortran_push_H,
  .pml_a                 = psc_push_fields_fortran_pml_a,
  .pml_b                 = psc_push_fields_fortran_pml_b,
};
