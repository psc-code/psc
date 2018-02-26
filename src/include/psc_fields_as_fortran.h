
#ifndef PSC_FIELDS_AS_FORTRAN_H
#define PSC_FIELDS_AS_FORTRAN_H

#include "psc_fields_fortran.h"

typedef PscMfieldsFortran     mfields_t;

#define psc_mfields_get_from          psc_mfields_fortran_get_from
#define psc_mfields_put_to            psc_mfields_fortran_put_to
#define FIELDS_TYPE                   "fortran"

#define PSC_FIELDS_AS_FORTRAN 1

#endif
