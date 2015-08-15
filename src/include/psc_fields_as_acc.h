
#ifndef PSC_FIELDS_AS_ACC_H
#define PSC_FIELDS_AS_ACC_H

#include "psc_fields_acc.h"

typedef struct psc_mfields mfields_t;
typedef struct psc_fields fields_t;
typedef fields_acc_real_t fields_real_t;
#define MPI_FIELDS_REAL MPI_FIELDS_ACC_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_ACC(pf, fldnr, jx,jy,jz)

#define FIELDS_TYPE                   "acc"

#define PSC_FIELDS_AS_ACC 1

#endif
