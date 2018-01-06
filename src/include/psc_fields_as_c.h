
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef fields_c_real_t fields_real_t;
typedef fields_c_t      fields_t;
typedef mfields_c_t     mfields_t;
#define fields_t_ctor            fields_c_t_ctor
#define fields_t_dtor            fields_c_t_dtor

#define MPI_FIELDS_REAL MPI_FIELDS_C_REAL

#define fields_t_set_nan              double_set_nan
#define FIELDS_TYPE                   "c"

#define PSC_FIELDS_AS_C 1

#endif
