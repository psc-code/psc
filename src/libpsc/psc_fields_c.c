
#include "psc.h"
#include "psc_fields_c.h"

#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void
__fields_c_alloc(fields_c_t *pf, int ib[3], int ie[3], int nr_comp,
		 fields_c_real_t *arr, bool with_array)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  if (with_array) {
    pf->flds = arr;
  } else {
#ifdef USE_CBE
    // The Cell processor translation can use the C fields with one modification:
    // the data needs to be 128 byte aligned (to speed off-loading to spes). This
    // change is roughly put in below.
    void *m;
    int ierr = posix_memalign(&m, 128, nr_comp * size * sizeof(*pf->flds));
    pf->flds =  m; 
    assert(ierr == 0);
#else
    pf->flds = calloc(nr_comp * size, sizeof(*pf->flds));
#endif
  }
  pf->with_array = with_array;
  pf->name = calloc(nr_comp, sizeof(*pf->name));
  for (int m = 0; m < nr_comp; m++) {
    pf->name[m] = NULL;
  }
}

void
fields_c_alloc(fields_c_t *pf, int ib[3], int ie[3], int nr_comp)
{
  __fields_c_alloc(pf, ib, ie, nr_comp, NULL, false);
}

void
fields_c_alloc_with_array(fields_c_t *pf, int ib[3], int ie[3], int nr_comp,
			  fields_c_real_t *arr)
{
  __fields_c_alloc(pf, ib, ie, nr_comp, arr, true);
}

void
fields_c_free(fields_c_t *pf)
{
  if (!pf->with_array) {
    free(pf->flds);
  }
  for (int m = 0; m < pf->nr_comp; m++) {
    free(pf->name[m]);
  }
  free(pf->name);
}

#if FIELDS_BASE == FIELDS_C

void
psc_mfields_c_get_from(mfields_c_t *flds, int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;
  *flds = *flds_base;
}

void
psc_mfields_c_put_to(mfields_c_t *flds, int mb, int me, void *_flds_base)
{
  flds->f = NULL;
}

#elif FIELDS_BASE == FIELDS_FORTRAN

void
psc_mfields_c_get_from(mfields_c_t *flds, int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_get", 1., 0, 0);
  }
  prof_start(pr);

  flds->domain = flds_base->domain;
  flds->nr_fields = flds_base->nr_fields;
  for (int d = 0; d < 3; d++) {
    flds->ibn[d] = flds_base->ibn[d];
  }
  psc_mfields_c_alloc(flds);

  psc_foreach_patch(ppsc, p) {
    fields_c_t *pf = &flds->f[p];
    fields_base_t *pf_base = &flds_base->f[p];
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_C(pf, m, jx,jy,jz) = F3_FORTRAN(pf_base, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }

  prof_stop(pr);
}

void
psc_mfields_c_put_to(mfields_c_t *flds, int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_c_put", 1., 0, 0);
  }
  prof_start(pr);

  psc_foreach_patch(ppsc, p) {
    fields_c_t *pf = &flds->f[p];
    fields_base_t *pf_base = &flds_base->f[p];
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_FORTRAN(pf_base, m, jx,jy,jz) = F3_C(pf, m, jx,jy,jz);
      }
    } foreach_3d_g_end;
  }

  psc_mfields_c_free(flds);

  prof_stop(pr);
}

#endif

void
fields_c_zero(fields_c_t *pf, int m)
{
  memset(&F3_C(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_c_real_t));
}

void
fields_c_set(fields_c_t *pf, int m, fields_c_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_C(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_c_copy_comp(fields_c_t *pto, int m_to, fields_c_t *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3_C(pto, m_to, jx, jy, jz) = F3_C(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

void
fields_c_axpy(fields_c_t *y, fields_c_real_t a, fields_c_t *x)
{
  assert(y->nr_comp == x->nr_comp);
  for (int m = 0; m < y->nr_comp; m++) {
    for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
      for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
	for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	  F3_C(y, m, jx, jy, jz) += a * F3_C(x, m, jx, jy, jz);
	}
      }
    }
  }
}

void
fields_c_scale(fields_c_t *pf, fields_c_real_t val)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
      for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
	for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	  F3_C(pf, m, jx, jy, jz) *= val;
	}
      }
    }
  }
}

void
psc_mfields_c_alloc(mfields_c_t *flds)
{
  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));
  for (int p = 0; p < flds->nr_patches; p++) {
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_c_alloc(&flds->f[p], ilg, ihg, flds->nr_fields);
  }
}

void
psc_mfields_c_free(mfields_c_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_c_free(&flds->f[p]);
  }
  free(flds->f);
  flds->f = NULL;
}

void
psc_mfields_c_zero(mfields_c_t *flds, int m)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_c_zero(&flds->f[p], m);
  }
}

void
psc_mfields_c_axpy(mfields_c_t *yf, fields_c_real_t alpha, mfields_c_t *xf)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_axpy(&yf->f[p], alpha, &xf->f[p]);
  }
}

void
psc_mfields_c_scale(mfields_c_t *yf, fields_c_real_t alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_scale(&yf->f[p], alpha);
  }
}

void
psc_mfields_c_set_comp(mfields_c_t *yf, int m, fields_c_real_t alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_set(&yf->f[p], m, alpha);
  }
}

void
psc_mfields_c_copy_comp(mfields_c_t *to, int mto, mfields_c_t *from, int mfrom)
{
  for (int p = 0; p < to->nr_patches; p++) {
    fields_c_copy_comp(&to->f[p], mto, &from->f[p], mfrom);
  }
}

