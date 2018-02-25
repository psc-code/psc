
#include "psc_push_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_E

static void
psc_push_fields_fortran_push_E(struct psc_push_fields *push,
			       struct psc_mfields *mflds_base)
{
  PscMfieldsFortran mf = mflds_base->get_as(JXI, HZ + 1);
  for (int p = 0; p < mf.nr_patches(); p++) {
    PIC_msa_e(psc_mfields_get_patch(mf.mflds(), p));
  }
  mf.put_as(mflds_base, EX, HZ + 1);
}

// ----------------------------------------------------------------------
// psc_push_fields_fortran_push_H

static void
psc_push_fields_fortran_push_H(struct psc_push_fields *push,
			       struct psc_mfields *mflds_base)
{
  PscMfieldsFortran mf = mflds_base->get_as(JXI, HZ + 1);
  for (int p = 0; p < mf.nr_patches(); p++) {
    PIC_msa_h(psc_mfields_get_patch(mf.mflds(), p));
  }
  mf.put_as(mflds_base, EX, HZ + 1);
}

// ======================================================================
// psc_push_fields: subclass "fortran"

struct psc_push_fields_ops psc_push_fields_fortran_ops = {
  .name                  = "fortran",
  .push_mflds_E          = psc_push_fields_fortran_push_mflds_E,
  .push_mflds_H          = psc_push_fields_fortran_push_mflds_H,
};
