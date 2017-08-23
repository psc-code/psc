
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
psc_mfields_fortran_copy_to_c(struct psc_mfields *mflds_fortran, struct psc_mfields *mflds_c,
			     int mb, int me)
{
  for (int p = 0; p < mflds_fortran->nr_patches; p++) {
    struct psc_fields *flds_fortran = psc_mfields_get_patch(mflds_fortran, p);
    struct psc_fields *flds_c = psc_mfields_get_patch(mflds_c, p);
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
}

static void
psc_mfields_fortran_copy_from_c(struct psc_mfields *mflds_fortran, struct psc_mfields *mflds_c,
			       int mb, int me)
{
  for (int p = 0; p < mflds_fortran->nr_patches; p++) {
    struct psc_fields *flds_fortran = psc_mfields_get_patch(mflds_fortran, p);
    struct psc_fields *flds_c = psc_mfields_get_patch(mflds_c, p);
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
}

static struct mrc_obj_method psc_mfields_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_mfields_fortran_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_mfields_fortran_copy_from_c),
  {}
};

#include "psc_fields_common.c"

