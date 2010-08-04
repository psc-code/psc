
#ifndef PSC_FIELD_FORTRAN_H
#define PSC_FIELD_FORTRAN_H

#include "psc.h"

typedef double fields_fortran_real_t;
#define MPI_FIELDS_FORTRAN_REAL MPI_DOUBLE

#if 1

#define F3_FORTRAN(fldnr, jx,jy,jz)		\
  (psc.f_fields[fldnr][FF3_OFF(jx,jy,jz)])

#else

#define F3_FORTRAN(fldnr, jx,jy,jz)					\
  (*({int off = FF3_OFF(jx,jy,jz);					\
      assert(off >= 0);							\
      assert(off < NR_FIELDS*psc.fld_size);				\
      &(psc.f_fields[fldnr][off]);					\
    }))

#endif

#endif
