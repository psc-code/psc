
#ifndef PSC_FIELDS_AS_CUDA2_H
#define PSC_FIELDS_AS_CUDA2_H

#include "psc_fields_cuda2.h"

typedef fields_cuda2_real_t fields_real_t;
#define MPI_FIELDS_REAL MPI_FIELDS_CUDA2_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_CUDA2(pf, fldnr, jx,jy,jz)

#define FIELDS_TYPE                   "cuda2"

#define PSC_FIELDS_AS_CUDA2 1

#endif
