
#include "psc_push_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_a_E

static void
psc_push_fields_fortran_push_a_E(struct psc_push_fields *push,
				 mfields_base_t *flds_base)
{
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
  
  PIC_msa_e(&flds.f[0]);
  
  fields_fortran_put(&flds, EX, HZ + 1, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_a_H

static void
psc_push_fields_fortran_push_a_H(struct psc_push_fields *push,
				 mfields_base_t *flds_base)
{
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
  
  PIC_msa_h(&flds.f[0]);
  
  fields_fortran_put(&flds, EX, HZ + 1, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_b_H

static void
psc_push_fields_fortran_push_b_H(struct psc_push_fields *push,
				 mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
  
  PIC_msb_h(&flds.f[0]);
    
  fields_fortran_put(&flds, EX, HZ + 1, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_b_E

static void
psc_push_fields_fortran_push_b_E(struct psc_push_fields *push,
				 mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
  
  PIC_msb_e(&flds.f[0]);
    
  fields_fortran_put(&flds, EX, HZ + 1, flds_base);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_pml_a

static void
psc_push_fields_fortran_pml_a(struct psc_push_fields *push,
			       mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
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
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_pml_b

static void
psc_push_fields_fortran_pml_b(struct psc_push_fields *push,
			      mfields_base_t *flds_base)
{
  assert(ppsc->nr_patches == 1);
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
}

// ======================================================================
// psc_push_fields: subclass "fortran"

struct psc_push_fields_ops psc_push_fields_fortran_ops = {
  .name                  = "fortran",
  .push_a_E              = psc_push_fields_fortran_push_a_E,
  .push_a_H              = psc_push_fields_fortran_push_a_H,
  .push_b_H              = psc_push_fields_fortran_push_b_H,
  .push_b_E              = psc_push_fields_fortran_push_b_E,
  .pml_a                 = psc_push_fields_fortran_pml_a,
  .pml_b                 = psc_push_fields_fortran_pml_b,
};
