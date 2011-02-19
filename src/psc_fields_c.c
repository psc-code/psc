
#include "psc.h"
#include "psc_fields_c.h"

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
    pf->flds = calloc(nr_comp * size, sizeof(*pf->flds));
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
  return __fields_c_alloc(pf, ib, ie, nr_comp, NULL, false);
}

void
fields_c_alloc_with_array(fields_c_t *pf, int ib[3], int ie[3], int nr_comp,
			  fields_c_real_t *arr)
{
  return __fields_c_alloc(pf, ib, ie, nr_comp, arr, true);
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
fields_c_get(fields_c_t *pf, int mb, int me)
{
  *pf = psc.pf;
}

void
fields_c_put(fields_c_t *pf, int mb, int me)
{
  pf->flds = NULL;
}

#else

void
fields_c_get(fields_c_t *pf, int mb, int me)
{
  fields_c_alloc(pf, psc.ilg, psc.ihg, NR_FIELDS);

  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_C(pf, m, jx+psc.ilo[0],jy+psc.ilo[1],jz+psc.ilo[2]) = F3_BASE(&psc.pf, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

void
fields_c_put(fields_c_t *pf, int mb, int me)
{
  foreach_patch(patch) {
    for (int m = mb; m < me; m++) {
      foreach_3d_g(patch, jx, jy, jz) {
	F3_BASE(&psc.pf, m, jx,jy,jz) = F3_C(pf, m, jx+psc.ilo[0],jy+psc.ilo[1],jz+psc.ilo[2]);
      }
    } foreach_3d_g_end;
  }

  fields_c_free(pf);
  pf->flds = NULL;
}

#endif

void
fields_c_zero(fields_c_t *pf, int m)
{
  memset(&F3_C(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_c_real_t));
}

void
fields_c_zero_all(fields_c_t *pf)
{
  memset(pf->flds, 0,
	 pf->nr_comp * pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_c_real_t));
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
fields_c_copy(fields_c_t *pf, int m_to, int m_from)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_C(pf, m_to, jx, jy, jz) = F3_C(pf, m_from, jx, jy, jz);
      }
    }
  }
}

void
fields_c_axpy_all(fields_c_t *y, fields_c_real_t a, fields_c_t *x)
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
fields_c_scale_all(fields_c_t *pf, fields_c_real_t val)
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

