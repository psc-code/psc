
#include "psc.h"
#include "psc_fields_cuda.h"
#include "psc_cuda.h"

#include <mrc_params.h>
#include <mrc_profile.h>

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
psc_mfields_copy_cf_to_cuda(mfields_cuda_t *flds_cuda, int mb, int me, mfields_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf_cuda = psc_mfields_get_patch_cuda(flds_cuda, p);
    fields_t *pf = psc_mfields_get_patch(flds, p);
    float *h_flds = calloc(12 * pf_cuda->im[0] * pf_cuda->im[1] * pf_cuda->im[2],
			   sizeof(*h_flds));

    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_CUDA(pf_cuda, m, jx,jy,jz) = F3(pf, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
    __fields_cuda_to_device(pf_cuda, h_flds, mb, me);
    free(h_flds);
  }
}

void
psc_mfields_copy_cf_from_cuda(mfields_cuda_t *flds_cuda, int mb, int me, mfields_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf_cuda = psc_mfields_get_patch_cuda(flds_cuda, p);
    fields_t *pf = psc_mfields_get_patch(flds, p);

    float *h_flds = calloc(12 * pf_cuda->im[0] * pf_cuda->im[1] * pf_cuda->im[2],
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
	F3(pf, m, jx,jy,jz) = F3_CUDA(pf_cuda, m, jx,jy,jz);
      }
    } foreach_3d_g_end;

    free(h_flds);
  }
}

static bool __gotten;

static struct psc_mfields *
_psc_mfields_cuda_get_c(struct psc_mfields *_flds_base, int mb, int me)
{
  mfields_cuda_t *flds_base = (mfields_cuda_t *) _flds_base;

  assert(!__gotten);
  __gotten = true;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_c_t *flds = psc_mfields_create(psc_comm(ppsc));
  psc_mfields_set_type(flds, "c");
  psc_mfields_set_domain(flds, flds_base->domain);
  psc_mfields_set_param_int(flds, "nr_fields", flds_base->nr_fields);
  psc_mfields_set_param_int3(flds, "ibn", ppsc->ibn);
  psc_mfields_setup(flds);
  psc_mfields_copy_cf_from_cuda(flds_base, mb, me, flds);

  prof_stop(pr);

  return (struct psc_mfields *) flds;
}

static void
_psc_mfields_cuda_put_c(struct psc_mfields *flds, struct psc_mfields *_flds_base, int mb, int me)
{
  assert(__gotten);
  __gotten = false;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_put", 1., 0, 0);
  }
  prof_start(pr);

  mfields_cuda_t *flds_base = (mfields_cuda_t *) _flds_base;
  psc_mfields_copy_cf_to_cuda(flds_base, mb, me, (mfields_c_t *) flds);
  psc_mfields_destroy(flds);

  prof_stop(pr);
}

struct psc_mfields *
_psc_mfields_c_get_cuda(struct psc_mfields *_flds_base, int mb, int me)
{
  mfields_c_t *flds_base = (mfields_c_t *) _flds_base;
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_cuda_t *flds = psc_mfields_create(psc_comm(ppsc));
  psc_mfields_set_type(flds, "cuda");
  psc_mfields_set_domain(flds, flds_base->domain);
  psc_mfields_set_param_int(flds, "nr_fields", flds_base->nr_fields);
  psc_mfields_set_param_int3(flds, "ibn", ppsc->ibn);
  psc_mfields_setup(flds);

  psc_mfields_copy_cf_to_cuda(flds, mb, me, flds_base);

  prof_stop(pr);

  return (struct psc_mfields *) flds;
}

void
_psc_mfields_c_put_cuda(struct psc_mfields *flds, struct psc_mfields *_flds_base, int mb, int me)
{
  mfields_c_t *flds_base = (mfields_c_t *) _flds_base;
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_put", 1., 0, 0);
  }
  prof_start(pr);

  psc_mfields_copy_cf_from_cuda((mfields_cuda_t *) flds, mb, me, flds_base);
  psc_mfields_destroy(flds);

  prof_stop(pr);
}

// ======================================================================
// psc_fields_cuda

static void
_psc_mfields_cuda_setup(mfields_cuda_t *flds)
{
  flds->data = calloc(ppsc->nr_patches, sizeof(fields_cuda_t));

  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = psc_mfields_get_patch_cuda(flds, p);
    struct psc_patch *patch = &ppsc->patch[p];
    for (int d = 0; d < 3; d++) {
      pf->ib[d] = -ppsc->ibn[d];
      pf->im[d] = patch->ldims[d] + 2 * ppsc->ibn[d];
    }

    pf->name = calloc(flds->nr_fields, sizeof(*pf->name));
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
    
    for (int m = 0; m < pf->nr_comp; m++) {
      free(pf->name[m]);
    }
    free(pf->name);
  }
  
  free(flds->data);
}

static struct psc_mfields *
_psc_mfields_cuda_get_cuda(struct psc_mfields *base, int mb, int me)
{
  return base;
}

static void
_psc_mfields_cuda_put_cuda(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me)
{
}

// ======================================================================
// psc_mfields: subclass "cuda"
  
struct psc_mfields_ops psc_mfields_cuda_ops = {
  .name                  = "cuda",
  .setup                 = _psc_mfields_cuda_setup,
  .destroy               = _psc_mfields_cuda_destroy,
  .get_c                 = _psc_mfields_cuda_get_c,
  .put_c                 = _psc_mfields_cuda_put_c,
  .get_cuda              = _psc_mfields_cuda_get_cuda,
  .put_cuda              = _psc_mfields_cuda_put_cuda,
};


