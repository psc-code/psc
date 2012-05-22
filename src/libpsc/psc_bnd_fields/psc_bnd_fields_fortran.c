
#include "psc_bnd_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fields_fortran_fill_ghosts_b_H

static void
psc_bnd_fields_fortran_fill_ghosts_b_H(struct psc_bnd_fields *bnd,
				       struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "fortran", JXI, HZ + 1);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_H_b", 1., 0, 0);
  }
  prof_start(pr);
  PIC_fill_ghosts_h_b(ppsc, flds->p, flds);
  prof_stop(pr);
  
  psc_fields_put_as(flds, flds_base, HX, HZ + 1);
}

// ======================================================================
// psc_bnd_fields: subclass "fortran"

struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops = {
  .name                  = "fortran",
  .fill_ghosts_b_H       = psc_bnd_fields_fortran_fill_ghosts_b_H,
};
