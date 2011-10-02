
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef mfields_c_t mfields_t;
typedef fields_c_t fields_t;
typedef fields_c_real_t fields_real_t;
#define MPI_FIELDS_REAL MPI_FIELDS_C_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_C(pf, fldnr, jx,jy,jz)

#define psc_mfields_get_from          psc_mfields_c_get_from
#define psc_mfields_put_to  	      psc_mfields_c_put_to
#define psc_mfields_get_patch         psc_mfields_c_get_patch_c
#define fields_size                   fields_c_size

#endif
