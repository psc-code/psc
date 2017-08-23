
#include "psc.h"
#include "psc_fields_as_single.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

#include "psc_fields_common.c"

// ======================================================================
// convert to c

static void
psc_fields_single_copy_from_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			      int mb, int me)
{
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

void
psc_fields_single_copy_to_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			    int mb, int me)
{
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

// ======================================================================

// ======================================================================
// psc_mfields: subclass "single"
  
struct psc_mfields_ops psc_mfields_single_ops = {
  .name                  = "single",
};

// ======================================================================
// psc_fields: subclass "single"
  
static struct mrc_obj_method psc_fields_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_single_copy_from_c),
  {}
};

struct psc_fields_ops psc_fields_single_ops = {
  .name                  = "single",
  .methods               = psc_fields_single_methods,
  .setup                 = psc_fields_single_setup,
  .destroy               = psc_fields_single_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_single_write,
  .read                  = psc_fields_single_read,
#endif
  .zero_comp             = psc_fields_single_zero_comp,
  .set_comp              = psc_fields_single_set_comp,
  .scale_comp            = psc_fields_single_scale_comp,
  .copy_comp             = psc_fields_single_copy_comp,
  .axpy_comp             = psc_fields_single_axpy_comp,
};

