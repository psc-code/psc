
#include "psc.h"
#include "psc_fields_fortran.h"

void
fields_fortran_alloc(fields_fortran_t *pf)
{
  f_real **fields = ALLOC_field();
  for (int i = 0; i < NR_FIELDS; i++) {
    pf->flds[i] = fields[i];
  }
}

void
fields_fortran_free(fields_fortran_t *pf)
{
  FREE_field();
  for (int i = 0; i < NR_FIELDS; i++) {
    pf->flds[i] = NULL;
  }
}

#if FIELDS_BASE == FIELDS_FORTRAN

void
fields_fortran_get(fields_fortran_t *pf, int mb, int me)
{
  fields_fortran_t *pf_base = &psc.pf;
  for (int i = 0; i < NR_FIELDS; i++) {
    pf->flds[i] = pf_base->flds[i];
  }
}

void
fields_fortran_put(fields_fortran_t *pf, int mb, int me)
{
}

#else

static fields_fortran_t __flds;
static int __gotten; // to check we're pairing get/put correctly

void
fields_fortran_get(fields_fortran_t *pf, int mb, int me)
{
  assert(!__gotten);
  __gotten = 1;

  if (!__flds.flds[0]) {
    fields_fortran_alloc(&__flds);
  }
  for (int m = 0; m < NR_FIELDS; m++) {
    pf->flds[m] = __flds.flds[m];
  }

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_FORTRAN(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }
}

void
fields_fortran_put(fields_fortran_t *pf, int mb, int me)
{
  assert(__gotten);
  __gotten = 0;

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_BASE(m, jx,jy,jz) = F3_FORTRAN(pf, m, jx,jy,jz);
	}
      }
    }
  }

  for (int m = 0; m < NR_FIELDS; m++) {
    pf->flds[0] = NULL;
  }
}

#endif

void
fields_fortran_zero(fields_fortran_t *pf, int m)
{
  memset(pf->flds[m], 0, psc.fld_size * sizeof(f_real));
}

void
fields_fortran_set(fields_fortran_t *pf, int m, fields_fortran_real_t val)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_FORTRAN(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_fortran_copy(fields_fortran_t *pf, int m_to, int m_from)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_FORTRAN(pf, m_to, jx, jy, jz) = F3_FORTRAN(pf, m_from, jx, jy, jz);
      }
    }
  }
}
