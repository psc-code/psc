
#include "psc.h"
#include "psc_fields_cuda.h"
#include "psc_cuda.h"

#include <mrc_params.h>

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

static void
psc_fields_cuda_setup(struct psc_fields *pf)
{
  __fields_cuda_alloc(pf);
}

static void
psc_fields_cuda_destroy(struct psc_fields *pf)
{
  __fields_cuda_free(pf);
}

#include "psc_fields_as_c.h"

static void
psc_fields_cuda_copy_from_c(struct psc_fields *flds_cuda, struct psc_fields *flds_c,
			    int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_CUDA(flds_cuda, m, jx,jy,jz) = F3(flds_c, m, jx,jy,jz);
	}
      }
    }
  }

  __fields_cuda_to_device(flds_cuda, h_flds, mb, me);

  free(h_flds);
}

void
psc_fields_cuda_copy_to_c(struct psc_fields *flds_cuda, struct psc_fields *flds_c,
			  int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  __fields_cuda_from_device(flds_cuda, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
#if 0
	  if (isnan(F3_CUDA(pf_cuda, m, jx,jy,jz))) {
	    printf("m %d j %d,%d,%d %g\n", m, jx,jy,jz,
		   F3_CUDA(pf_cuda, m, jx,jy,jz));
	    assert(0);
	  }
#endif
	  F3_C(flds_c, m, jx,jy,jz) = F3_CUDA(flds_cuda, m, jx,jy,jz);
	}
      }
    }
  }

  free(h_flds);
}

// ======================================================================
// psc_mfields: subclass "cuda"
  
struct psc_mfields_ops psc_mfields_cuda_ops = {
  .name                  = "cuda",
};

// ======================================================================
// psc_fields: subclass "cuda"
  
static struct mrc_obj_method psc_fields_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_cuda_copy_from_c),
  {}
};

struct psc_fields_ops psc_fields_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_fields_cuda),
  .methods               = psc_fields_cuda_methods,
  .setup                 = psc_fields_cuda_setup,
  .destroy               = psc_fields_cuda_destroy,
};

