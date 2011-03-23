
#include "psc_bnd_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fields_fortran_fill_ghosts_b_H

static void
psc_bnd_fields_fortran_fill_ghosts_b_H(struct psc_bnd_fields *bnd,
				       mfields_base_t *flds_base)
{
  assert(psc.nr_patches == 1);
  mfields_fortran_t flds;
  fields_fortran_get(&flds, JXI, HZ + 1, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_H_b", 1., 0, 0);
  }
  prof_start(pr);
  PIC_fill_ghosts_h_b(&flds.f[0]);
  prof_stop(pr);
  
  fields_fortran_put(&flds, HX, HZ + 1, flds_base);
}

// ======================================================================
// psc_bnd_fields: subclass "fortran"

struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops = {
  .name                  = "fortran",
  .fill_ghosts_b_H       = psc_bnd_fields_fortran_fill_ghosts_b_H,
};
