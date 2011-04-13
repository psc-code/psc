
#ifndef PSC_FIELDS_AS_FORTRAN_H
#define PSC_FIELDS_AS_FORTRAN_H

#include "psc_fields_fortran.h"

typedef fields_fortran_t fields_t;
typedef mfields_fortran_t mfields_t;

#define F3(fldnr, jx,jy,jz) F3_FORTRAN(pf, fldnr, jx,jy,jz)

#define fields_get          fields_fortran_get
#define fields_put  	    fields_fortran_put
#define fields_zero         fields_fortran_zero

#endif
