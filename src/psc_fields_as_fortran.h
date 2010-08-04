
#ifndef PSC_FIELDS_AS_FORTRAN_H
#define PSC_FIELDS_AS_FORTRAN_H

#include "psc_fields_fortran.h"

typedef psc_fields_fortran_t psc_fields_t;

#define F3(fldnr, jx,jy,jz) F3_FORTRAN(pf, fldnr, jx,jy,jz)

#define psc_fields_get          psc_fields_fortran_get
#define psc_fields_put  	psc_fields_fortran_put
#define psc_fields_zero         psc_fields_fortran_zero

#endif
