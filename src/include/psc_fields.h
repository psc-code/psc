
#ifndef PSC_FIELDS_H
#define PSC_FIELDS_H

// ----------------------------------------------------------------------
// mfields type

// This type is replicated for each actual fields type, however,
// the interface and implementation is always identical, hence 
// created automatically for the variants using macros

#define DECLARE_MFIELDS_METHODS(type)					\
  									\
typedef struct psc_mfields_##type {				        \
  struct mrc_obj obj;							\
  fields_##type##_t *f;							\
  int nr_patches;							\
  list_t entry;								\
  struct mrc_domain *domain;						\
  int nr_fields;							\
  int ibn[3];								\
} mfields_##type##_t;							\
									\
MRC_CLASS_DECLARE(psc_mfields_##type, struct psc_mfields_##type);	\
									\
void psc_mfields_##type##_set_domain(mfields_##type##_t *flds,	        \
				     struct mrc_domain *domain);	\
void psc_mfields_##type##_get_from(mfields_##type##_t *pf, int mb, int me, void *flds_base); \
void psc_mfields_##type##_put_to(mfields_##type##_t *pf, int mb, int me, void *flds_base); \
									\
/* FIXME, should be per mrc_domain or sth, really */			\
extern list_t mfields_##type##_list;					\


#include "psc_fields_fortran.h"
DECLARE_MFIELDS_METHODS(fortran)

#include "psc_fields_c.h"
DECLARE_MFIELDS_METHODS(c)

#include "psc_fields_sse2.h"
DECLARE_MFIELDS_METHODS(sse2)

// ----------------------------------------------------------------------
// base fields type

#if FIELDS_BASE == FIELDS_FORTRAN

typedef fields_fortran_t fields_base_t;
typedef fields_fortran_real_t fields_base_real_t;
typedef mfields_fortran_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_FORTRAN_REAL

#define fields_base_alloc            fields_fortran_alloc
#define fields_base_alloc_with_array fields_fortran_alloc_with_array
#define fields_base_free             fields_fortran_free
#define fields_base_zero             fields_fortran_zero
#define fields_base_zero_all         fields_fortran_zero_all
#define fields_base_set              fields_fortran_set
#define fields_base_copy             fields_fortran_copy
#define fields_base_axpy_all         fields_fortran_axpy_all
#define fields_base_scale_all        fields_fortran_scale_all
#define fields_base_size             fields_fortran_size
#define mfields_base_list            mfields_fortran_list
#define psc_mfields_base_create         psc_mfields_fortran_create
#define psc_mfields_base_set_name       psc_mfields_fortran_set_name
#define psc_mfields_base_set_param_int  psc_mfields_fortran_set_param_int
#define psc_mfields_base_set_param_int3 psc_mfields_fortran_set_param_int3
#define psc_mfields_base_setup          psc_mfields_fortran_setup
#define psc_mfields_base_destroy        psc_mfields_fortran_destroy
#define psc_mfields_base_write          psc_mfields_fortran_write
#define psc_mfields_base_read           psc_mfields_fortran_read
#define psc_mfields_base_set_domain     psc_mfields_fortran_set_domain


#define F3_BASE(pf, m, jx,jy,jz)     F3_FORTRAN(pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_C

#include "psc_fields_c.h"

typedef fields_c_t fields_base_t;
typedef fields_c_real_t fields_base_real_t;
typedef mfields_c_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_C_REAL

#define fields_base_alloc            fields_c_alloc
#define fields_base_alloc_with_array fields_c_alloc_with_array
#define fields_base_free             fields_c_free
#define fields_base_zero             fields_c_zero
#define fields_base_zero_all         fields_c_zero_all
#define fields_base_set              fields_c_set
#define fields_base_copy             fields_c_copy
#define fields_base_axpy_all         fields_c_axpy_all
#define fields_base_scale_all        fields_c_scale_all
#define fields_base_size             fields_c_size
#define mfields_base_list            mfields_c_list
#define psc_mfields_base_create         psc_mfields_c_create
#define psc_mfields_base_set_name       psc_mfields_c_set_name
#define psc_mfields_base_set_param_int  psc_mfields_c_set_param_int
#define psc_mfields_base_set_param_int3 psc_mfields_c_set_param_int3
#define psc_mfields_base_setup          psc_mfields_c_setup
#define psc_mfields_base_destroy        psc_mfields_c_destroy
#define psc_mfields_base_write          psc_mfields_c_write
#define psc_mfields_base_read           psc_mfields_c_read
#define psc_mfields_base_set_domain     psc_mfields_c_set_domain

#define F3_BASE(pf, m, jx,jy,jz)     F3_C(pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_SSE2

#include "psc_fields_sse2.h"

typedef fields_sse2_t fields_base_t;
typedef fields_sse2_real_t fields_base_real_t;
typedef mfields_sse2_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL MPI_FIELDS_SSE2_REAL

#define fields_base_alloc fields_sse2_alloc
#define fields_base_free  fields_sse2_free
#define fields_base_zero  fields_sse2_zero
#define fields_base_set   fields_sse2_set
#define fields_base_copy  fields_sse2_copy
#define mfields_base_list            mfields_sse2_list
#define psc_mfields_base_create      psc_mfields_sse2_create
#define psc_mfields_base_setup       psc_mfields_sse2_setup
#define psc_mfields_base_destroy     psc_mfields_sse2_destroy
#define psc_mfields_base_set_domain  psc_mfields_sse2_set_domain

#define F3_BASE(pf, m, jx,jy,jz) F3_SSE2(pf, m, jx,jy,jz)

#else
#error unknown FIELDS_BASE
#endif

#endif


