
#ifndef PSC_FIELD_SSE2_H
#define PSC_FIELD_SSE2_H

#include "psc.h"
#include "simd_sse2.h"

typedef sse2_real fields_sse2_real_t;

typedef struct {
  fields_sse2_real_t *flds;
} fields_sse2_t;

// FIXME, this needs to be looked into for efficiency
#define F3_OFF_SSE2(fldnr, jx,jy,jz)					\
  ((((((fldnr)								\
       *(psc.patch[0].ldims[2] + 2*psc.ibn[0]) + ((jz)+psc.ibn[2]))	\
      *(psc.patch[0].ldims[1] + 2*psc.ibn[1]) + ((jy)+psc.ibn[1]))	\
     *(psc.patch[0].ldims[0] + 2*psc.ibn[0]) + ((jx)+psc.ibn[0]))))

#if 1

#define F3_SSE2(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_SSE2(fldnr, jx,jy,jz)])

#else
//out of range debugging
#define F3_SSE2(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_SSE2(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * psc.fld_size);				\
      &((pf)->flds[off]);						\
    }))

#endif

void fields_sse2_alloc(fields_sse2_t *pf);
void fields_sse2_free(fields_sse2_t *pf);
void fields_sse2_get(fields_sse2_t *pf, int mb, int me, void *);
void fields_sse2_put(fields_sse2_t *pf, int mb, int me, void *);
void fields_sse2_zero(fields_sse2_t *pf, int m);

#endif
