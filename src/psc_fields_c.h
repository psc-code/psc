
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc.h"

typedef double fields_c_real_t;
#define MPI_FIELDS_C_REAL MPI_DOUBLE

typedef struct {
  fields_c_real_t *flds;
} fields_c_t;

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

void fields_c_alloc(fields_c_t *pf);
void fields_c_free(fields_c_t *pf);
void fields_c_get(fields_c_t *pf, int mb, int me);
void fields_c_put(fields_c_t *pf, int mb, int me);
void fields_c_zero(fields_c_t *pf, int m);
void fields_c_set(fields_c_t *pf, int m, fields_c_real_t val);
void fields_c_copy(fields_c_t *pf, int m_to, int m_from);

#endif
