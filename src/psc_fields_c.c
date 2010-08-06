
#include "psc.h"
#include "psc_fields_c.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// FIXME, we're allocating all fields even if the PML ones aren't
// going to be used

void
fields_c_alloc(fields_c_t *pf)
{
  pf->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));
}

void
fields_c_free(fields_c_t *pf)
{
  free(pf->flds);
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

#elif FIELDS_BASE == FIELDS_FORTRAN

static fields_c_real_t *__flds;
static int __gotten;

void
fields_c_get(fields_c_t *pf, int mb, int me)
{
  assert(!__gotten);
  __gotten = 1;

  if (!__flds) {
    __flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));
  }
  pf->flds = __flds;

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_C(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }
}

void
fields_c_put(fields_c_t *pf, int mb, int me)
{
  assert(__gotten);
  __gotten = 0;

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_BASE(m, jx,jy,jz) = F3_C(pf, m, jx,jy,jz);
	}
      }
    }
  }

  pf->flds = NULL;
}

#endif

void
fields_c_zero(fields_c_t *pf, int m)
{
  memset(&F3_C(pf, m, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(fields_c_real_t));
}

void
fields_c_set(fields_c_t *pf, int m, fields_c_real_t val)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_C(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_c_copy(fields_c_t *pf, int m_to, int m_from)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_C(pf, m_to, jx, jy, jz) = F3_C(pf, m_from, jx, jy, jz);
      }
    }
  }
}
