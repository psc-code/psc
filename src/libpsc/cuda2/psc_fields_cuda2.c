
#include "psc_fields_cuda2.h"

// for conversions
#include "psc_fields_single.h"

#include <stdlib.h>

// ======================================================================
// psc_fields "cuda2"

// ----------------------------------------------------------------------
// psc_fields_cuda2_setup

static void
psc_fields_cuda2_setup(struct psc_fields *pf)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    size *= pf->im[d];
  }
  pf->data = calloc(pf->nr_comp * size, sizeof(fields_cuda2_real_t));
}

// ----------------------------------------------------------------------
// psc_fields_cuda2_destroy

static void
psc_fields_cuda2_destroy(struct psc_fields *pf)
{
  free(pf->data);
}

// ----------------------------------------------------------------------
// convert from/to "single"

static void
psc_fields_cuda2_copy_from_single(struct psc_fields *flds_cuda2, struct psc_fields *flds_single,
				 int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda2->ib[2]; jz < flds_cuda2->ib[2] + flds_cuda2->im[2]; jz++) {
      for (int jy = flds_cuda2->ib[1]; jy < flds_cuda2->ib[1] + flds_cuda2->im[1]; jy++) {
	for (int jx = flds_cuda2->ib[0]; jx < flds_cuda2->ib[0] + flds_cuda2->im[0]; jx++) {
	  F3_CUDA2(flds_cuda2, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	}
      }
    }
  }

  //  __fields_cuda2_to_device(flds_cuda2, h_flds, mb, me);
}

static void
psc_fields_cuda2_copy_to_single(struct psc_fields *flds_cuda2, struct psc_fields *flds_single,
			       int mb, int me)
{
  //  __fields_cuda2_from_device(flds_cuda2, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda2->ib[2]; jz < flds_cuda2->ib[2] + flds_cuda2->im[2]; jz++) {
      for (int jy = flds_cuda2->ib[1]; jy < flds_cuda2->ib[1] + flds_cuda2->im[1]; jy++) {
	for (int jx = flds_cuda2->ib[0]; jx < flds_cuda2->ib[0] + flds_cuda2->im[0]; jx++) {
	  F3_S(flds_single, m, jx,jy,jz) = F3_CUDA2(flds_cuda2, m, jx,jy,jz);
	}
      }
    }
  }
}

// ======================================================================
// psc_fields: subclass "cuda2"
  
static struct mrc_obj_method psc_fields_cuda2_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_fields_cuda2_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_fields_cuda2_copy_from_single),
  {}
};

struct psc_fields_ops psc_fields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_fields_cuda2),
  .methods               = psc_fields_cuda2_methods,
  .setup                 = psc_fields_cuda2_setup,
  .destroy               = psc_fields_cuda2_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_cuda2_write,
#endif
  .axpy_comp             = psc_fields_cuda2_axpy_comp,
  .zero_comp             = psc_fields_cuda2_zero_comp,
#endif
};

// ======================================================================
// psc_mfields "cuda2"

// ======================================================================
// psc_mfields: subclass "cuda2"
  
struct psc_mfields_ops psc_mfields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_mfields_cuda2),
#if 0
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_mfields_cuda2_read,
#endif
#endif
};
