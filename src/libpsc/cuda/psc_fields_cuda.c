
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

static void
fields_cuda_alloc(mfields_cuda_t *flds, int nr_comp)
{
  flds->f = calloc(ppsc->nr_patches, sizeof(*flds->f));

  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = &flds->f[p];
    struct psc_patch *patch = &ppsc->patch[p];
    for (int d = 0; d < 3; d++) {
      pf->ib[d] = -ppsc->ibn[d];
      pf->im[d] = patch->ldims[d] + 2 * ppsc->ibn[d];
    }

    pf->nr_comp = nr_comp;
    pf->name = calloc(nr_comp, sizeof(*pf->name));

    __fields_cuda_alloc(pf);
  }
}

static void
fields_cuda_free(mfields_cuda_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf = &flds->f[p];
    __fields_cuda_free(pf);
    
    for (int m = 0; m < pf->nr_comp; m++) {
      free(pf->name[m]);
    }
    free(pf->name);
  }
  
  free(flds->f);
  flds->f = NULL;
}

#include "psc_fields_as_c.h"

void
psc_mfields_copy_cf_to_cuda(mfields_cuda_t *flds_cuda, int mb, int me, mfields_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_cuda_t *pf_cuda = &flds_cuda->f[p];
    fields_t *pf = &flds->f[p];
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
    fields_cuda_t *pf_cuda = &flds_cuda->f[p];
    fields_t *pf = &flds->f[p];

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

#if FIELDS_BASE == FIELDS_CUDA

mfields_cuda_t *
psc_mfields_cuda_get_from(int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;
  return flds_base;
}

void
psc_mfields_cuda_put_to(mfields_cuda_t *flds, int mb, int me, void *_flds_base)
{
}

static bool __gotten;

mfields_c_t *
psc_mfields_c_get_from(int mb, int me, void *_flds_base)
{
  mfields_cuda_t *flds_base = _flds_base;

  assert(!__gotten);
  __gotten = true;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_c_t *flds = psc_mfields_c_create(psc_comm(ppsc));
  psc_mfields_c_set_domain(flds, flds_base->domain);
  psc_mfields_c_set_param_int(flds, "nr_fields", flds_base->nr_fields);
  psc_mfields_c_set_param_int3(flds, "ibn", ppsc->ibn);
  psc_mfields_c_setup(flds);
  psc_mfields_copy_cf_from_cuda(flds_base, mb, me, flds);

  prof_stop(pr);

  return flds;
}

void
psc_mfields_c_put_to(mfields_c_t *flds, int mb, int me, void *_flds_base)
{
  assert(__gotten);
  __gotten = false;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_put", 1., 0, 0);
  }
  prof_start(pr);

  mfields_cuda_t *flds_base = _flds_base;
<<<<<<< HEAD
  psc_mfields_copy_cf_to_cuda(flds_base, mb, me, flds);
  psc_mfields_c_free(flds);
  free(flds);
=======
  psc_mfields_copy_cf_from_cuda(flds_base, mb, me, flds);
  psc_mfields_c_destroy(flds);
>>>>>>> 630cd2b... psc/mfields: get rid of psc_mfields_alloc, use psc_mfields_create

  prof_stop(pr);
}

mfields_fortran_t *
psc_mfields_fortran_get_from(int mb, int me, void *_flds_base)
{
  assert(0);
}

void
psc_mfields_fortran_put_to(mfields_fortran_t *flds, int mb, int me, void *_flds_base)
{
  assert(0);
}

#elif FIELDS_BASE == FIELDS_C

mfields_cuda_t *
psc_mfields_cuda_get_from(int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_cuda_t *flds = calloc(1, sizeof(*flds));
  fields_cuda_alloc(flds, 12);
  psc_mfields_copy_cf_to_cuda(flds, mb, me, flds_base);

  prof_stop(pr);

  return flds;
}

void
psc_mfields_cuda_put_to(mfields_cuda_t *flds, int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;
  static int pr;
  if (!pr) {
    pr = prof_register("fields_cuda_put", 1., 0, 0);
  }
  prof_start(pr);

  psc_mfields_copy_cf_from_cuda(flds, mb, me, flds_base);
  fields_cuda_free(flds);
  free(flds);

  prof_stop(pr);
}

#endif

// ======================================================================
// psc_mfields_cuda_list

LIST_HEAD(psc_mfields_cuda_list);

void
psc_mfields_cuda_list_add(mfields_cuda_t **flds_p)
{
  mfields_cuda_list_entry_t *p = malloc(sizeof(*p));
  p->flds_p = flds_p;
  list_add_tail(&p->entry, &psc_mfields_cuda_list);
}

void
psc_mfields_cuda_list_del(mfields_cuda_t **flds_p)
{
  mfields_cuda_list_entry_t *p;
  __list_for_each_entry(p, &psc_mfields_cuda_list, entry, mfields_cuda_list_entry_t) {
    if (p->flds_p == flds_p) {
      list_del(&p->entry);
      free(p);
      return;
    }
  }
  assert(0);
}

// ======================================================================
// psc_fields_cuda

LIST_HEAD(mfields_cuda_list);

#define VAR(x) (void *)offsetof(struct psc_mfields_cuda, x)
static struct param psc_mfields_cuda_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  {},
};
#undef VAR

void
psc_mfields_cuda_set_domain(mfields_cuda_t *flds, struct mrc_domain *domain)
{
  flds->domain = domain;
}

static void
_psc_mfields_cuda_setup(mfields_cuda_t *flds)
{
  fields_cuda_alloc(flds, flds->nr_fields);
}

static void
_psc_mfields_cuda_destroy(mfields_cuda_t *flds)
{
  fields_cuda_free(flds);
}

void
psc_mfields_cuda_axpy(mfields_cuda_t *yf, fields_cuda_real_t alpha, mfields_cuda_t *xf)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    assert(0);
    //    fields_cuda_axpy(&yf->f[p], alpha, &xf->f[p]);
  }
}

void
psc_mfields_cuda_scale(mfields_cuda_t *yf, fields_cuda_real_t alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    assert(0);
    //    fields_cuda_scale(&yf->f[p], alpha);
  }
}

void
psc_mfields_cuda_set_comp(mfields_cuda_t *yf, int m, fields_cuda_real_t alpha)
{
  assert(0);
}

void
psc_mfields_cuda_copy_comp(mfields_cuda_t *to, int mto, mfields_cuda_t *from, int mfrom)
{
  assert(0);
}

struct mrc_class_psc_mfields_cuda mrc_class_psc_mfields_cuda = {
  .name             = "psc_mfields_cuda",
  .size             = sizeof(struct psc_mfields_cuda),
  .param_descr      = psc_mfields_cuda_descr,
  .setup            = _psc_mfields_cuda_setup,
  .destroy          = _psc_mfields_cuda_destroy,
};

