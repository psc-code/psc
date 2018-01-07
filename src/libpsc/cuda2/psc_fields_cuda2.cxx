
#include "psc_fields_cuda2.h"
#include "psc_cuda2.h"

// for conversions
#include "psc_fields_single.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// psc_fields "cuda2"

// ----------------------------------------------------------------------
// psc_fields_cuda2_zero_comp

static void
psc_fields_cuda2_zero_comp(struct psc_fields *pf, int m)
{
  memset(&F3_CUDA2(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_cuda2_real_t));
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

// ----------------------------------------------------------------------
// convert from/to "cuda"

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_CUDA(pf, fldnr, jx,jy,jz)				\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_CUDA(pf, fldnr, jx,jy,jz)		\
  (h_flds[F3_OFF_CUDA(pf, fldnr, jx,jy,jz)])

#else

#define F3_CUDA(pf, fldnr, jx,jy,jz)				\
  (*({int off = F3_OFF_CUDA(pf, fldnr, jx,jy,jz);			\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(h_flds[off]);						\
    }))

#endif

#ifdef USE_CUDA

EXTERN_C void __fields_cuda_to_device(struct psc_fields *pf, real *h_flds, int mb, int me);
EXTERN_C void __fields_cuda_from_device(struct psc_fields *pf, real *h_flds, int mb, int me);

static void
psc_fields_cuda2_copy_from_cuda(struct psc_fields *flds_cuda2, struct psc_fields *flds_cuda,
				int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  __fields_cuda_from_device(flds_cuda, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_CUDA2(flds_cuda2, m, jx,jy,jz) = F3_CUDA(flds_cuda, m, jx,jy,jz);
	}
      }
    }
  }

  free(h_flds);
}

static void
psc_fields_cuda2_copy_to_cuda(struct psc_fields *flds_cuda2, struct psc_fields *flds_cuda,
			      int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_CUDA(flds_cuda, m, jx,jy,jz) = F3_CUDA2(flds_cuda2, m, jx,jy,jz);
	}
      }
    }
  }

  __fields_cuda_to_device(flds_cuda, h_flds, mb, me);

  free(h_flds);
}

#endif

// ======================================================================
// psc_fields: subclass "cuda2"
  
static struct mrc_obj_method psc_fields_cuda2_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_fields_cuda2_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_fields_cuda2_copy_from_single),
#ifdef USE_CUDA
  MRC_OBJ_METHOD("copy_to_cuda"    , psc_fields_cuda2_copy_to_cuda),
  MRC_OBJ_METHOD("copy_from_cuda"  , psc_fields_cuda2_copy_from_cuda),
#endif
  {}
};

struct psc_fields_ops psc_fields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_fields_cuda2),
  .methods               = psc_fields_cuda2_methods,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_cuda2_write,
#endif
  .axpy_comp             = psc_fields_cuda2_axpy_comp,
#endif
  .zero_comp             = psc_fields_cuda2_zero_comp,
};

// ======================================================================
// psc_mfields "cuda2"

// ----------------------------------------------------------------------
// psc_mfields_cuda2_setup

static void
psc_mfields_cuda2_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

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

  sub->h_flds = calloc(mflds->nr_fields * total_size, sizeof(*sub->h_flds));

  fields_cuda2_real_t *h_flds = sub->h_flds;
  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    unsigned int size = flds->im[0] * flds->im[1] * flds->im[2];
    flds->data = h_flds;
    h_flds += mflds->nr_fields * size;
  }

  sub->d_flds = cuda_calloc(mflds->nr_fields * total_size, sizeof(*sub->d_flds));
}

// ----------------------------------------------------------------------
// psc_mfields_cuda2_destroy

static void
psc_mfields_cuda2_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

    flds->data = NULL;
  }

  free(sub->h_flds);
  cuda_free(sub->d_flds);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda2_copy_to_device

void
psc_mfields_cuda2_copy_to_device(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

  size_t total_size = sub->im[0] * sub->im[1] * sub->im[2] * mflds->nr_patches * mflds->nr_fields;
  cuda_memcpy_device_from_host(sub->d_flds, sub->h_flds, total_size * sizeof(*sub->d_flds));
}

// ----------------------------------------------------------------------
// psc_mfields_cuda2_copy_to_host

void
psc_mfields_cuda2_copy_to_host(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

  size_t total_size = sub->im[0] * sub->im[1] * sub->im[2] * mflds->nr_patches * mflds->nr_fields;
  cuda_memcpy_host_from_device(sub->h_flds, sub->d_flds, total_size * sizeof(*sub->d_flds));
}

// ======================================================================
// psc_mfields: subclass "cuda2"
  
struct psc_mfields_ops psc_mfields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_mfields_cuda2),
  .setup                 = psc_mfields_cuda2_setup,
  .destroy               = psc_mfields_cuda2_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_mfields_cuda2_read,
#endif
#endif
};
