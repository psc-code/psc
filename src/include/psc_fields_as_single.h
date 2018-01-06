
#ifndef PSC_FIELDS_AS_SINGLE_H
#define PSC_FIELDS_AS_SINGLE_H

#include "psc_fields_single.h"

typedef fields_single_real_t fields_real_t;
typedef fields_single_t      fields_t;
typedef mfields_single_t     mfields_t;
#define fields_t_mflds           fields_single_t_mflds
#define fields_t_size            fields_single_t_size
#define fields_t_zero_range      fields_single_t_zero_range

#define MPI_FIELDS_REAL MPI_FIELDS_SINGLE_REAL

#define fields_t_set_nan              float_set_nan
#define FIELDS_TYPE                   "single"

#define PSC_FIELDS_AS_SINGLE 1

#endif
