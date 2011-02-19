
#include "psc.h"
#include "psc_fields_fortran.h"

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
    static bool ALLOC_field_called;
    if (!ALLOC_field_called && nr_comp == NR_FIELDS &&
	ib[0] == psc.ilg[0] && ib[1] == psc.ilg[1] && ib[2] == psc.ilg[2] &&
	ie[0] == psc.ihg[0] && ie[1] == psc.ihg[1] && ie[2] == psc.ihg[2]) {
      ALLOC_field_called = true;
      f_real **fields = ALLOC_field();
      for (int i = 0; i < NR_FIELDS; i++) {
	pf->flds[i] = fields[i];
      }
      pf->fortran_alloc = true;
    } else {
      for (int i = 0; i < nr_comp; i++) {
	pf->flds[i] = calloc(sizeof(*pf->flds[i]), size);
      }
      pf->fortran_alloc = false;
    }
  }
}

void
fields_fortran_alloc(fields_fortran_t *pf, int ib[3], int ie[3], int nr_comp)
{
  return __fields_fortran_alloc(pf, ib, ie, nr_comp, NULL, false);
}

void
fields_fortran_alloc_with_array(fields_fortran_t *pf, int ib[3], int ie[3],
				int nr_comp, fields_fortran_real_t *arr)
{
  return __fields_fortran_alloc(pf, ib, ie, nr_comp, arr, true);
}


void
fields_fortran_free(fields_fortran_t *pf)
{
  if (!pf->with_array) {
    if (pf->fortran_alloc) {
      FREE_field();
    } else {
      for (int i = 0; i < pf->nr_comp; i++) {
	free(pf->flds[i]);
      }
    }
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
fields_fortran_get(fields_fortran_t *pf, int mb, int me)
{
  fields_fortran_t *pf_base = &psc.pf;
  *pf = *pf_base;
}

void
fields_fortran_get_from(fields_fortran_t *pf, int mb, int me,
			void *_pf_base, int mb_base)
{
  fields_fortran_get(pf, mb, me);

  fields_base_t *pf_base = _pf_base;
  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_FORTRAN(pf, m, jx,jy,jz) = F3_BASE(pf_base, m - mb + mb_base, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

void
fields_fortran_put(fields_fortran_t *pf, int mb, int me)
{
}

void
fields_fortran_put_to(fields_fortran_t *pf, int mb, int me,
		      void *_pf_base, int mb_base)
{
  fields_base_t *pf_base = _pf_base;
  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_BASE(pf_base, m - mb + mb_base, jx,jy,jz) = F3_FORTRAN(pf, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
  
  fields_fortran_put(pf, mb, me);
}

#else

static fields_fortran_t __flds;
static int __gotten; // to check we're pairing get/put correctly

void
fields_fortran_get_from(fields_fortran_t *pf, int mb, int me,
			void *_pf_base, int mb_base)
{
  fields_base_t *pf_base = _pf_base;
  assert(!__gotten);
  __gotten = 1;

  if (!__flds.flds) {
    struct psc_patch *patch = &psc.patch[0];
    fields_fortran_alloc(&__flds, patch->ilg, patch->ihg, NR_FIELDS);
  }
  *pf = __flds;

  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_FORTRAN(pf, m, jx,jy,jz) = F3_BASE(pf_base, m - mb + mb_base, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

void
fields_fortran_get(fields_fortran_t *pf, int mb, int me)
{
  fields_fortran_get_from(pf, mb, me, &psc.pf, mb);
}

void
fields_fortran_put_to(fields_fortran_t *pf, int mb, int me,
		      void *_pf_base, int mb_base)
{
  fields_base_t *pf_base = _pf_base;
  assert(__gotten);
  __gotten = 0;

  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_BASE(pf_base, m - mb + mb_base, jx,jy,jz) = F3_FORTRAN(pf, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

void
fields_fortran_put(fields_fortran_t *pf, int mb, int me)
{
  return fields_fortran_put_to(pf, mb, me, &psc.pf, mb);
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
  foreach_patch(patch) {
    foreach_3d_g(patch, jx, jy, jz) {
      F3_FORTRAN(pf, m, jx,jy,jz) = val;
    } foreach_3d_g_end;
  }
}

void
fields_fortran_copy(fields_fortran_t *pf, int m_to, int m_from)
{
  foreach_patch(patch) {
    foreach_3d_g(patch, jx, jy, jz) {
      F3_FORTRAN(pf, m_to, jx,jy,jz) = F3_FORTRAN(pf, m_from, jx,jy,jz);
    } foreach_3d_g_end;
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

