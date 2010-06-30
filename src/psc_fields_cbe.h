
#ifndef PSC_FIELD_CBE_H
#define PSC_FIELD_CBE_H

#include "psc.h"
#include "simd_cbe.h"

#if CBE_DOUBLE
typedef double fields_cbe_real_t;
#define MPI_FIELD_CBE_REAL MPI_DOUBLE

#else 
typedef float fields_cbe_real_t;
#define MPI_FIELD_CBE_REAL MPI_FLOAT
#endif

typedef struct { 
  fields_cbe_real_t *flds;
} fields_cbe_t;

#define F3_OFF_CBE(fldnr, jx,jy,jz)					\
  ((((((jz)-psc.ilg[2]))						\
     *psc.img[1] + ((jy)-psc.ilg[1]))					\
    *psc.img[0] + ((jx)-psc.ilg[0]))					\
   *NR_FIELDS + fldnr)

#if 0

#define F3_CBE(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_CBE(fldnr, jx,jy,jz)])

#else
//out of range debugging
#define F3_CBE(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_CBE(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * psc.fld_size);				\
      &((pf)->flds[off]);						\
    }))

#endif

void fields_cbe_alloc(fields_cbe_t *pf);
void fields_cbe_free(fields_cbe_t *pf);
void fields_cbe_get(fields_cbe_t *pf, int mb, int me);
void fields_cbe_put(fields_cbe_t *pf, int mb, int mf);
void fields_cbe_zero(fields_cbe_t *pf, int m);
void fields_cbe_set(fields_cbe_t *pf, int m, fields_cbe_real_t val);
void fields_cbe_copy(fields_cbe_t *pf, int m_to, int m_from);

#endif
