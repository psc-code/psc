
#include "psc.h"
#include "psc_fields_fortran.h"
#include "psc_glue.h"

#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

void
__fields_fortran_alloc(fields_fortran_t *pf, int ib[3], int ie[3], int nr_comp,
		       fields_fortran_real_t *arr, bool with_array)
{
  pf->flds = calloc(nr_comp, sizeof(*pf->flds));
  pf->name = calloc(nr_comp, sizeof(*pf->name));
  for (int m = 0; m < nr_comp; m++) {
    pf->name[m] = NULL;
  }

  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->with_array = with_array;

  if (with_array) {
    assert(nr_comp == 1);
    for (int i = 0; i < nr_comp; i++) {
      pf->flds[i] = arr;
    }
  } else {
    pf->flds[0] = calloc(size * nr_comp, sizeof(*pf->flds[0]));
    for (int i = 1; i < nr_comp; i++) {
      pf->flds[i] = pf->flds[0] + i * size;
    }
  }
}

void
fields_fortran_alloc(fields_fortran_t *pf, int ib[3], int ie[3], int nr_comp)
{
  __fields_fortran_alloc(pf, ib, ie, nr_comp, NULL, false);
}

void
fields_fortran_alloc_with_array(fields_fortran_t *pf, int ib[3], int ie[3],
				int nr_comp, fields_fortran_real_t *arr)
{
  __fields_fortran_alloc(pf, ib, ie, nr_comp, arr, true);
}


void
fields_fortran_free(fields_fortran_t *pf)
{
  if (!pf->with_array) {
    free(pf->flds[0]);
  }
  for (int i = 0; i < pf->nr_comp; i++) {
    pf->flds[i] = NULL;
  }
  for (int m = 0; m < pf->nr_comp; m++) {
    free(pf->name[m]);
  }
  free(pf->name);
  free(pf->flds);
}

#if FIELDS_BASE == FIELDS_FORTRAN

void
psc_mfields_fortran_get_from(mfields_fortran_t *flds, int mb, int me, void *_flds_base)
{
  mfields_base_t *flds_base = _flds_base;
  *flds = *flds_base;
}

void
psc_mfields_fortran_put_to(mfields_fortran_t *flds, int mb, int me, void *_flds_base)
{
}

#elif FIELDS_BASE == FIELDS_C

void
psc_mfields_fortran_get_from(mfields_fortran_t *flds, int mb, int me, void *_flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fields_fortran_get", 1., 0, 0);
  }
  prof_start(pr);

  mfields_base_t *flds_base = _flds_base;
  flds->f = calloc(ppsc->nr_patches, sizeof(*flds->f));
  psc_foreach_patch(ppsc, p) {
    fields_fortran_t *pf = &flds->f[p];
    struct psc_patch *patch = &ppsc->patch[p];
    int ilg[3] = { -ppsc->ibn[0], -ppsc->ibn[1], -ppsc->ibn[2] };
    int ihg[3] = { patch->ldims[0] + ppsc->ibn[0],
		   patch->ldims[1] + ppsc->ibn[1],
		   patch->ldims[2] + ppsc->ibn[2] };
    fields_fortran_alloc(pf, ilg, ihg, NR_FIELDS);

    fields_base_t *pf_base = &flds_base->f[p];
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_FORTRAN(pf, m, jx,jy,jz) = F3_C(pf_base, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }

  prof_stop(pr);
}

void
psc_mfields_fortran_put_to(mfields_fortran_t *flds, int mb, int me, void *_flds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fields_fortran_put", 1., 0, 0);
  }
  prof_start(pr);

  mfields_base_t *flds_base = _flds_base;
  psc_foreach_patch(ppsc, p) {
    fields_fortran_t *pf = &flds->f[p];
    fields_base_t *pf_base = &flds_base->f[p];
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_C(pf_base, m, jx,jy,jz) = F3_FORTRAN(pf, m, jx,jy,jz);
      }
    } foreach_3d_g_end;

    fields_fortran_free(pf);
  }
  
  free(flds->f);
  flds->f = NULL;

  prof_stop(pr);
}

#endif

void
fields_fortran_zero(fields_fortran_t *pf, int m)
{
  memset(pf->flds[m], 0, 
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_fortran_real_t));
}

void
fields_fortran_zero_all(fields_fortran_t *pf)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    fields_fortran_zero(pf, m);
  }
}

void
fields_fortran_set(fields_fortran_t *pf, int m, fields_fortran_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_FORTRAN(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_fortran_copy(fields_fortran_t *pf, int m_to, int m_from)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_FORTRAN(pf, m_to, jx,jy,jz) = F3_FORTRAN(pf, m_from, jx,jy,jz);
      }
    }
  }
}

void
fields_fortran_axpy_all(fields_fortran_t *y, fields_fortran_real_t a,
			fields_fortran_t *x)
{
  assert(y->nr_comp == x->nr_comp);
  for (int m = 0; m < y->nr_comp; m++) {
    for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
      for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
	for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	  F3_FORTRAN(y, m, jx, jy, jz) += a * F3_FORTRAN(x, m, jx, jy, jz);
	}
      }
    }
  }
}

void
fields_fortran_scale_all(fields_fortran_t *pf, fields_fortran_real_t val)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
      for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
	for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	  F3_FORTRAN(pf, m, jx, jy, jz) *= val;
	}
      }
    }
  }
}

