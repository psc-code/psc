
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc.h"

typedef double psc_fields_c_real_t;

typedef struct {
  psc_fields_c_real_t *flds;
} psc_fields_c_t;

#define F3_OFF_C(fldnr, jx,jy,jz)					\
  (((((fldnr								\
       *psc.img[2] + ((jz)-psc.ilg[2]))					\
      *psc.img[1] + ((jy)-psc.ilg[1]))					\
     *psc.img[0] + ((jx)-psc.ilg[0]))))

#if 1

#define F3_C(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_C(fldnr, jx,jy,jz)])

#else

#define F3_C(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_C(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS*psc.fld_size);				\
      &((pf)->flds[off]);						\
    }))

#endif

void psc_fields_c_get(psc_fields_c_t *pf, int mb, int me);
void psc_fields_c_put(psc_fields_c_t *pf, int mb, int me);

#endif
