
#ifndef PSC_FIELDS_AS_FORTRAN_H
#define PSC_FIELDS_AS_FORTRAN_H

#include "psc_fields_fortran.h"

typedef fields_fortran_real_t fields_real_t;

#define F3(fldnr, jx,jy,jz) F3_FORTRAN(pf, fldnr, jx,jy,jz)

#define psc_mfields_get_from          psc_mfields_fortran_get_from
#define psc_mfields_put_to            psc_mfields_fortran_put_to

#define PSC_FIELDS_AS_FORTRAN 1

#endif
