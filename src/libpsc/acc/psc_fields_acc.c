
#include "psc_fields_acc.h"
#include "psc_acc.h"

// for conversions
#include "psc_fields_single.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// psc_fields "acc"

// ----------------------------------------------------------------------
// psc_fields_acc_zero_comp

static void
psc_fields_acc_zero_comp(struct psc_fields *pf, int m)
{
  memset(&F3_ACC(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_acc_real_t));
}

// ----------------------------------------------------------------------
// convert from/to "single"

static void
psc_fields_acc_copy_from_single(struct psc_fields *flds_acc, struct psc_fields *flds_single,
				 int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_acc->ib[2]; jz < flds_acc->ib[2] + flds_acc->im[2]; jz++) {
      for (int jy = flds_acc->ib[1]; jy < flds_acc->ib[1] + flds_acc->im[1]; jy++) {
	for (int jx = flds_acc->ib[0]; jx < flds_acc->ib[0] + flds_acc->im[0]; jx++) {
	  F3_ACC(flds_acc, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	}
      }
    }
  }

  //  __fields_acc_to_device(flds_acc, h_flds, mb, me);
}

static void
psc_fields_acc_copy_to_single(struct psc_fields *flds_acc, struct psc_fields *flds_single,
			       int mb, int me)
{
  //  __fields_acc_from_device(flds_acc, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_acc->ib[2]; jz < flds_acc->ib[2] + flds_acc->im[2]; jz++) {
      for (int jy = flds_acc->ib[1]; jy < flds_acc->ib[1] + flds_acc->im[1]; jy++) {
	for (int jx = flds_acc->ib[0]; jx < flds_acc->ib[0] + flds_acc->im[0]; jx++) {
	  F3_S(flds_single, m, jx,jy,jz) = F3_ACC(flds_acc, m, jx,jy,jz);
	}
      }
    }
  }
}

// ======================================================================
// psc_fields: subclass "acc"
  
static struct mrc_obj_method psc_fields_acc_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_fields_acc_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_fields_acc_copy_from_single),
  {}
};

struct psc_fields_ops psc_fields_acc_ops = {
  .name                  = "acc",
  .size                  = sizeof(struct psc_fields_acc),
  .methods               = psc_fields_acc_methods,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_acc_write,
#endif
  .axpy_comp             = psc_fields_acc_axpy_comp,
#endif
  .zero_comp             = psc_fields_acc_zero_comp,
};

// ======================================================================
// psc_mfields "acc"

// ----------------------------------------------------------------------
// psc_mfields_acc_setup

static void
psc_mfields_acc_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_acc *sub = psc_mfields_acc(mflds);

  psc_mfields_setup_super(mflds);

  unsigned int total_size = 0;
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (p == 0) {
      for (int d = 0; d < 3; d++) {
	sub->im[d] = flds->im[d];
	sub->ib[d] = flds->ib[d];
      }
    } else {
      for (int d = 0; d < 3; d++) {
	assert(sub->im[d] == flds->im[d]);
	assert(sub->ib[d] == flds->ib[d]);
      }
    }

    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    total_size += size;
  }

  sub->flds = calloc(mflds->nr_fields * total_size, sizeof(*sub->flds));

  fields_acc_real_t *data = sub->flds;
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    flds->data = data;
    data += mflds->nr_fields * size;
  }
}

// ----------------------------------------------------------------------
// psc_mfields_acc_destroy

static void
psc_mfields_acc_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_acc *sub = psc_mfields_acc(mflds);

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

    flds->data = NULL;
  }

  free(sub->flds);
}

// ======================================================================
// psc_mfields: subclass "acc"
  
struct psc_mfields_ops psc_mfields_acc_ops = {
  .name                  = "acc",
  .size                  = sizeof(struct psc_mfields_acc),
  .setup                 = psc_mfields_acc_setup,
  .destroy               = psc_mfields_acc_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_mfields_acc_read,
#endif
#endif
};
