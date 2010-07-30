
#ifndef PSC_FIELD_SSE2_H
#define PSC_FIELD_SSE2_H

#include "psc.h"

#include <xmmintrin.h>

////////
/// Toggle to switch precision 
///
/// 1 is double precision, 0 is single. 

#define SSE2_DOUBLE 1

// ======================================================================

#if SSE2_DOUBLE

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 2

/// SSE2 floating point type
typedef double sse2_real;

#else

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 4 

/// SSE2 floating point type
typedef float sse2_real;

#endif

// ======================================================================

typedef sse2_real psc_fields_sse2_real_t;

#define MPI_F3_SSE2_REAL MPI_SSE2_REAL

struct psc_fields_sse2 {
  sse2_real *flds;
};

typedef struct psc_fields_sse2 psc_fields_sse2_t;

#define F3_OFF_SSE2(fldnr, jx,jy,jz)					\
  (((((fldnr								\
       *psc.img[2] + ((jz)-psc.ilg[2]))					\
      *psc.img[1] + ((jy)-psc.ilg[1]))					\
     *psc.img[0] + ((jx)-psc.ilg[0]))))

#if 1

#define F3_SSE2(pf, fldnr, jx,jy,jz)		\
  (pf->flds[F3_OFF_SSE2(fldnr, jx,jy,jz)])

#else
//out of range debugging
#define F3_SSE2(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_SSE2(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * psc.fld_size);				\
      &(pf->flds[off]);							\
    }))

#endif

//void psc_fields_sse2_alloc(struct psc_fields_c *pf);
//void psc_fields_sse2_free(struct psc_fields_c *pf);
void psc_fields_sse2_get(psc_fields_sse2_t *pf, int mb, int me);
void psc_fields_sse2_put(psc_fields_sse2_t *pf, int mb, int me);
//void psc_fields_sse2_zero(struct psc_fields_c *pf, int m);

#endif
