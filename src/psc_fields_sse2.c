
#include "psc.h"
#include "psc_fields_sse2.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if FIELDS_BASE == FIELDS_SSE2

void
fields_sse2_alloc(fields_sse2_t *pf)
{
  pf->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));
}

void
fields_sse2_free(fields_sse2_t *pf)
{
  free(pf->flds);
}

void
fields_sse2_get(fields_sse2_t *pf, int mb, int me)
{
  *pf = psc.pf;
}

void
fields_sse2_put(fields_sse2_t *pf, int mb, int me)
{
  pf->flds = NULL;
}

#else

static bool __gotten; // to check we're pairing get/put correctly

/// Copy fields from base data structure to an SSE2 friendly format.
void
fields_sse2_get(fields_sse2_t *pf, int mb, int me)
{
  assert(!__gotten);
  __gotten = true;

  pf->flds = _mm_malloc(NR_FIELDS*psc.fld_size*sizeof(sse2_real), 16);
  
  int *ilg = psc.ilg;
  for(int m = mb; m < me; m++){
    for(int n = 0; n < psc.fld_size; n++){
      //preserve Fortran ordering for now
      pf->flds[m * psc.fld_size + n] =
	(sse2_real) ((&F3_BASE(m, ilg[0],ilg[1],ilg[2]))[n]);
    }
  }
}

/// Copy fields from SSE2 data structures into base structures.
void
fields_sse2_put(fields_sse2_t *pf, int mb, int me)
{
  assert(__gotten);
  __gotten = false;

  int *ilg = psc.ilg;
  for(int m = mb; m < me; m++){
    for(int n = 0; n < psc.fld_size; n++){
      ((&F3_BASE(m, ilg[0],ilg[1],ilg[2]))[n]) = 
	pf->flds[m * psc.fld_size + n];
    }
  }
  _mm_free(pf->flds);
}

#endif

void
fields_sse2_zero(fields_sse2_t *pf, int m)
{
  memset(&F3_SSE2(pf, m, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(fields_sse2_real_t));
}

void
fields_sse2_set(fields_sse2_t *pf, int m, fields_sse2_real_t val)
{
  foreach_3d_g(jx, jy, jz) {
    F3_SSE2(pf, m, jx, jy, jz) = val;
  } foreach_3d_g_end;
}

void
fields_sse2_copy(fields_sse2_t *pf, int m_to, int m_from)
{
  foreach_3d_g(jx, jy, jz) {
    F3_SSE2(pf, m_to, jx, jy, jz) = F3_SSE2(pf, m_from, jx, jy, jz);
  } foreach_3d_g_end;
}
