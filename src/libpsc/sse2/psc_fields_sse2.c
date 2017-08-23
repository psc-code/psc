
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
psc_mfields_sse2_get_from(fields_sse2_t *pf, int mb, int me)
{
  *pf = psc.pf;
}

void
psc_mfields_sse2_put_to(fields_sse2_t *pf, int mb, int me)
{
  pf->flds = NULL;
}

#elif FIELDS_BASE == FIELDS_C

static bool __gotten; // to check we're pairing get/put correctly

/// Copy fields from base data structure to an SSE2 friendly format.
void
psc_mfields_sse2_get_from(fields_sse2_t *pf, int mb, int me, void *_flds_base)
{
  assert(!__gotten);
  __gotten = true;

  struct psc_patch *patch = &psc.patch[0];
  struct psc_mfields *flds_base = _flds_base;
  fields_base_t *pf_base = &flds_base->f[0];
  int sz = 1;
  for (int d = 0; d < 3; d++) {
    sz *= patch->ldims[d] + 2 * psc.ibn[d];
  }
  pf->flds = _mm_malloc(NR_FIELDS*sz*sizeof(sse2_real), 16);
  
  int *ibn = psc.ibn;
  for(int m = mb; m < me; m++){
    for(int n = 0; n < sz; n++){
      //preserve Fortran ordering for now
      pf->flds[m * sz + n] =
	(sse2_real) ((&F3_C(pf_base, m, -ibn[0],-ibn[1],-ibn[2]))[n]);
    }
  }
}

/// Copy fields from SSE2 data structures into base structures.
void
psc_mfields_sse2_put_to(fields_sse2_t *pf, int mb, int me, void *_flds_base)
{
  assert(__gotten);
  __gotten = false;

  struct psc_patch *patch = &psc.patch[0];
  struct psc_mfields *flds_base = _flds_base;
  fields_base_t *pf_base = &flds_base->f[0];
  int sz = 1;
  for (int d = 0; d < 3; d++) {
    sz *= patch->ldims[d] + 2 * psc.ibn[d];
  }
  int *ibn = psc.ibn;
  for(int m = mb; m < me; m++){
    for(int n = 0; n < sz; n++){
      ((&F3_C(pf_base, m, -ibn[0],-ibn[1],-ibn[2]))[n]) = 
	pf->flds[m * sz + n];
    }
  }
  _mm_free(pf->flds);
}

#endif

void
fields_sse2_zero(fields_sse2_t *pf, int m)
{
  struct psc_patch *patch = &psc.patch[0];
  int sz = 1;
  for (int d = 0; d < 3; d++) {
    sz *= patch->ldims[d] + 2 * psc.ibn[d];
  }
  memset(&F3_SSE2(pf, m, -psc.ibn[0], -psc.ibn[1], -psc.ibn[2]), 0,
	 sz * sizeof(fields_sse2_real_t));
}

void
fields_sse2_set(fields_sse2_t *pf, int m, fields_sse2_real_t val)
{
  foreach_patch(patch) {
    foreach_3d_g(patch, jx, jy, jz) {
      F3_SSE2(pf, m, jx,jy,jz) = val;
    } foreach_3d_g_end;
  }
}

void
fields_sse2_copy(fields_sse2_t *pf, int m_to, int m_from)
{
  foreach_patch(patch) {
    foreach_3d_g(patch, jx, jy, jz) {
      F3_SSE2(pf, m_to, jx,jy,jz) = F3_SSE2(pf, m_from, jx,jy,jz);
    } foreach_3d_g_end;
  }
}
