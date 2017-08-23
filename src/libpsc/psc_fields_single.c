
#include "psc.h"
#include "psc_fields_as_single.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

// ======================================================================
// convert to c

static void
psc_mfields_single_copy_from_c(struct psc_mfields *mflds_single, struct psc_mfields *mflds_c,
			       int mb, int me)
{
  for (int p = 0; p < mflds_single->nr_patches; p++) {
    struct psc_fields *flds_single = psc_mfields_get_patch(mflds_single, p);
    struct psc_fields *flds_c = psc_mfields_get_patch(mflds_c, p);
    for (int m = mb; m < me; m++) {
      for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
	for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	  for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	    F3_S(flds_single, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	  }
	}
      }
    }
  }
}

void
psc_mfields_single_copy_to_c(struct psc_mfields *mflds_single, struct psc_mfields *mflds_c,
			     int mb, int me)
{
  for (int p = 0; p < mflds_single->nr_patches; p++) {
    struct psc_fields *flds_single = psc_mfields_get_patch(mflds_single, p);
    struct psc_fields *flds_c = psc_mfields_get_patch(mflds_c, p);
    for (int m = mb; m < me; m++) {
      for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
	for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	  for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	    F3_C(flds_c, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	  }
	}
      }
    }
  }
}

static struct mrc_obj_method psc_mfields_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_mfields_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_mfields_single_copy_from_c),
  {}
};

#include "psc_fields_common.c"

