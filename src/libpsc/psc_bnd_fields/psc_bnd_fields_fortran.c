
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
  mfields_fortran_t *flds = psc_mfields_fortran_get_from(JXI, HZ + 1, flds_base);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_H_b", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    PIC_fill_ghosts_h_b(ppsc, p, &flds->f[p]);
  }
  prof_stop(pr);
  
  psc_mfields_fortran_put_to(flds, HX, HZ + 1, flds_base);
}

// ======================================================================
// psc_bnd_fields: subclass "fortran"

struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops = {
  .name                  = "fortran",
  .fill_ghosts_b_H       = psc_bnd_fields_fortran_fill_ghosts_b_H,
};
