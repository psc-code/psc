
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

#include "psc_fields_as_c.h"

void
psc_mfields_cuda_copy_from_c(mfields_cuda_t *flds_cuda, mfields_c_t *flds_c, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf_cuda = psc_mfields_get_patch_cuda(flds_cuda, p);
    fields_t *pf_c = psc_mfields_get_patch_c(flds_c, p);
    float *h_flds = calloc(flds_cuda->nr_fields * pf_cuda->im[0] * pf_cuda->im[1] * pf_cuda->im[2],
			   sizeof(*h_flds));

    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_CUDA(pf_cuda, m, jx,jy,jz) = F3(pf_c, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
    __fields_cuda_to_device(pf_cuda, h_flds, mb, me);
    free(h_flds);
  }
}

void
psc_mfields_cuda_copy_to_c(mfields_cuda_t *flds_cuda, mfields_c_t *flds_c, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf_cuda = psc_mfields_get_patch_cuda(flds_cuda, p);
    struct psc_fields *pf_c = psc_mfields_get_patch(flds_c, p);

    float *h_flds = calloc(flds_cuda->nr_fields * pf_cuda->im[0] * pf_cuda->im[1] * pf_cuda->im[2],
			   sizeof(*h_flds));
    __fields_cuda_from_device(pf_cuda, h_flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
#if 0
	if (isnan(F3_CUDA(pf_cuda, m, jx,jy,jz))) {
	  printf("m %d j %d,%d,%d %g\n", m, jx,jy,jz,
		 F3_CUDA(pf_cuda, m, jx,jy,jz));
	  assert(0);
	}
#endif
	F3_C(pf_c, m, jx,jy,jz) = F3_CUDA(pf_cuda, m, jx,jy,jz);
      } foreach_3d_g_end;
    }

    free(h_flds);
  }
}

// ======================================================================
// psc_fields_cuda

EXTERN_C void cuda_init(int rank);

static void
_psc_mfields_cuda_setup(mfields_cuda_t *flds)
{
  psc_mfields_setup_super(flds);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cuda_init(rank);

  flds->flds = calloc(ppsc->nr_patches, sizeof(*flds->flds));

  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = calloc(1, sizeof(*pf));
    flds->flds[p] = (struct psc_fields *) pf;
    struct psc_patch *patch = &ppsc->patch[p];
    for (int d = 0; d < 3; d++) {
      pf->ib[d] = -ppsc->ibn[d];
      pf->im[d] = patch->ldims[d] + 2 * ppsc->ibn[d];
    }
    pf->nr_comp = flds->nr_fields;

    __fields_cuda_alloc(pf);
  }
}

static void
_psc_mfields_cuda_destroy(mfields_cuda_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);
    __fields_cuda_free(pf);
    free(pf);
  }
  
  free(flds->flds);
}

// ======================================================================
// psc_mfields: subclass "cuda"
  
struct psc_mfields_ops psc_mfields_cuda_ops = {
  .name                  = "cuda",
  .setup                 = _psc_mfields_cuda_setup,
  .destroy               = _psc_mfields_cuda_destroy,
  .copy_to_c             = psc_mfields_cuda_copy_to_c,
  .copy_from_c           = psc_mfields_cuda_copy_from_c,
};



