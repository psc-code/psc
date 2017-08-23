
#include "psc.h"
#include "psc_fields_as_fortran.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

#include "psc_fields_inc.h"

// ======================================================================
// convert to/from "c"

static void
psc_fields_fortran_copy_to_c(struct psc_fields *flds_fortran, struct psc_fields *flds_c,
			     int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_c->ib[2]; jz < flds_c->ib[2] + flds_c->im[2]; jz++) {
      for (int jy = flds_c->ib[1]; jy < flds_c->ib[1] + flds_c->im[1]; jy++) {
	for (int jx = flds_c->ib[0]; jx < flds_c->ib[0] + flds_c->im[0]; jx++) {
	  F3_C(flds_c, m, jx,jy,jz) = F3_FORTRAN(flds_fortran, m, jx,jy,jz);
	}
      }
    }
  }
}

static void
psc_fields_fortran_copy_from_c(struct psc_fields *flds_fortran, struct psc_fields *flds_c,
			       int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_c->ib[2]; jz < flds_c->ib[2] + flds_c->im[2]; jz++) {
      for (int jy = flds_c->ib[1]; jy < flds_c->ib[1] + flds_c->im[1]; jy++) {
	for (int jx = flds_c->ib[0]; jx < flds_c->ib[0] + flds_c->im[0]; jx++) {
	  F3_FORTRAN(flds_fortran, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	}
      }
    }
  }
}

static struct mrc_obj_method psc_fields_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_fortran_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_fortran_copy_from_c),
  {}
};

#include "psc_fields_common.c"

// ======================================================================
// psc_mfields: subclass "fortran"
  
struct psc_mfields_ops psc_mfields_fortran_ops = {
  .name                  = "fortran",
};

// ======================================================================
// psc_fields: subclass "fortran"
  
struct psc_fields_ops psc_fields_fortran_ops = {
  .name                  = "fortran",
  .methods               = psc_fields_fortran_methods,
  .setup                 = psc_fields_fortran_setup,
  .destroy               = psc_fields_fortran_destroy,
  .zero_comp             = psc_fields_fortran_zero_comp,
  .set_comp              = psc_fields_fortran_set_comp,
  .scale_comp            = psc_fields_fortran_scale_comp,
  .copy_comp             = psc_fields_fortran_copy_comp,
  .axpy_comp             = psc_fields_fortran_axpy_comp,
};

