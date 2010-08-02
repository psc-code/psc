
#ifndef PSC_FIELD_FORTRAN_H
#define PSC_FIELD_FORTRAN_H

#include "psc.h"

typedef double fields_fortran_real_t;
#define MPI_FIELDS_FORTRAN_REAL MPI_DOUBLE

typedef struct {
  fields_fortran_real_t *flds[NR_FIELDS];
} psc_fields_fortran_t;

#if 1

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)                \
  ((pf)->flds[fldnr][FF3_OFF(jx,jy,jz)])

#else

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)					\
  (*({int off = FF3_OFF(jx,jy,jz);					\
      assert(off >= 0);							\
      assert(off < NR_FIELDS*psc.fld_size);				\
      &((pf)->flds[fldnr][off]);					\
    }))

#endif

#endif
