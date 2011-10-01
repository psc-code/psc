
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
void psc_mfields_##type##_axpy(mfields_##type##_t *yf, fields_##type##_real_t alpha, \
			       mfields_##type##_t *xf);			\
void psc_mfields_##type##_scale(mfields_##type##_t *yf, fields_##type##_real_t alpha); \
void psc_mfields_##type##_zero(mfields_##type##_t *flds, int m);	\
void psc_mfields_##type##_set_comp(mfields_##type##_t *flds, int m, fields_##type##_real_t alpha); \
void psc_mfields_##type##_copy_comp(mfields_##type##_t *to, int mto,	\
				    mfields_##type##_t *from, int mfrom);	\
void psc_mfields_##type##_alloc(mfields_##type##_t *flds_comp);	\
void psc_mfields_##type##_free(mfields_##type##_t *flds);		\
void psc_mfields_##type##_list_add(mfields_##type##_t **flds_p);	\
void psc_mfields_##type##_list_del(mfields_##type##_t **flds_p);	\
									\
typedef struct {							\
  mfields_##type##_t **flds_p;						\
  list_t entry;								\
} mfields_##type##_list_entry_t;					\
									\
/* FIXME, should be per mrc_domain or sth, really */			\
extern list_t mfields_##type##_list;					\
extern list_t psc_mfields_##type##_list;				\


#include "psc_fields_fortran.h"
DECLARE_MFIELDS_METHODS(fortran)

#include "psc_fields_c.h"
DECLARE_MFIELDS_METHODS(c)

#ifdef USE_SSE2
#include "psc_fields_sse2.h"
DECLARE_MFIELDS_METHODS(sse2)
#endif

#ifdef USE_CUDA
#include "psc_fields_cuda.h"
DECLARE_MFIELDS_METHODS(cuda)
#endif

// ----------------------------------------------------------------------
// base fields type

#if FIELDS_BASE == FIELDS_FORTRAN

typedef fields_fortran_t fields_base_t;
typedef fields_fortran_real_t fields_base_real_t;
typedef mfields_fortran_t mfields_base_t;
typedef mfields_fortran_list_entry_t mfields_base_list_entry_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_FORTRAN_REAL

#define mfields_base_list            mfields_fortran_list
#define psc_mfields_base_list           psc_mfields_fortran_list
#define psc_mfields_base_list_add       psc_mfields_fortran_list_add
#define psc_mfields_base_list_del       psc_mfields_fortran_list_del
#define psc_mfields_base_create         psc_mfields_fortran_create
#define psc_mfields_base_set_name       psc_mfields_fortran_set_name
#define psc_mfields_base_set_param_int  psc_mfields_fortran_set_param_int
#define psc_mfields_base_set_param_int3 psc_mfields_fortran_set_param_int3
#define psc_mfields_base_setup          psc_mfields_fortran_setup
#define psc_mfields_base_destroy        psc_mfields_fortran_destroy
#define psc_mfields_base_write          psc_mfields_fortran_write
#define psc_mfields_base_read           psc_mfields_fortran_read
#define psc_mfields_base_set_domain     psc_mfields_fortran_set_domain
#define psc_mfields_base_axpy           psc_mfields_fortran_axpy
#define psc_mfields_base_scale          psc_mfields_fortran_scale
#define psc_mfields_base_set_comp       psc_mfields_fortran_set_comp
#define psc_mfields_base_copy_comp      psc_mfields_fortran_copy_comp

#elif FIELDS_BASE == FIELDS_C

#include "psc_fields_c.h"

typedef fields_c_t fields_base_t;
typedef fields_c_real_t fields_base_real_t;
typedef mfields_c_t mfields_base_t;
typedef mfields_c_list_entry_t mfields_base_list_entry_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_C_REAL

#define mfields_base_list            mfields_c_list
#define psc_mfields_base_list           psc_mfields_c_list
#define psc_mfields_base_list_add       psc_mfields_c_list_add
#define psc_mfields_base_list_del       psc_mfields_c_list_del
#define psc_mfields_base_create         psc_mfields_c_create
#define psc_mfields_base_set_name       psc_mfields_c_set_name
#define psc_mfields_base_set_param_int  psc_mfields_c_set_param_int
#define psc_mfields_base_set_param_int3 psc_mfields_c_set_param_int3
#define psc_mfields_base_setup          psc_mfields_c_setup
#define psc_mfields_base_destroy        psc_mfields_c_destroy
#define psc_mfields_base_write          psc_mfields_c_write
#define psc_mfields_base_read           psc_mfields_c_read
#define psc_mfields_base_set_domain     psc_mfields_c_set_domain
#define psc_mfields_base_axpy           psc_mfields_c_axpy
#define psc_mfields_base_scale          psc_mfields_c_scale
#define psc_mfields_base_set_comp       psc_mfields_c_set_comp
#define psc_mfields_base_copy_comp      psc_mfields_c_copy_comp

#elif FIELDS_BASE == FIELDS_SSE2

#include "psc_fields_sse2.h"

typedef fields_sse2_t fields_base_t;
typedef fields_sse2_real_t fields_base_real_t;
typedef mfields_sse2_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL MPI_FIELDS_SSE2_REAL

#define mfields_base_list            mfields_sse2_list
#define psc_mfields_base_create      psc_mfields_sse2_create
#define psc_mfields_base_setup       psc_mfields_sse2_setup
#define psc_mfields_base_destroy     psc_mfields_sse2_destroy
#define psc_mfields_base_set_domain  psc_mfields_sse2_set_domain

#elif FIELDS_BASE == FIELDS_CUDA

#include "psc_fields_cuda.h"

typedef fields_cuda_t fields_base_t;
typedef fields_cuda_real_t fields_base_real_t;
typedef mfields_cuda_t mfields_base_t;
typedef mfields_cuda_list_entry_t mfields_base_list_entry_t;

#define MPI_FIELDS_BASE_REAL MPI_FIELDS_CUDA_REAL

#define mfields_base_list            mfields_cuda_list
#define psc_mfields_base_list           psc_mfields_cuda_list
#define psc_mfields_base_list_add       psc_mfields_cuda_list_add
#define psc_mfields_base_list_del       psc_mfields_cuda_list_del
#define psc_mfields_base_create         psc_mfields_cuda_create
#define psc_mfields_base_set_name       psc_mfields_cuda_set_name
#define psc_mfields_base_set_param_int  psc_mfields_cuda_set_param_int
#define psc_mfields_base_set_param_int3 psc_mfields_cuda_set_param_int3
#define psc_mfields_base_setup          psc_mfields_cuda_setup
#define psc_mfields_base_destroy        psc_mfields_cuda_destroy
#define psc_mfields_base_write          psc_mfields_cuda_write
#define psc_mfields_base_read           psc_mfields_cuda_read
#define psc_mfields_base_set_domain     psc_mfields_cuda_set_domain
#define psc_mfields_base_axpy           psc_mfields_cuda_axpy
#define psc_mfields_base_scale          psc_mfields_cuda_scale
#define psc_mfields_base_set_comp       psc_mfields_cuda_set_comp
#define psc_mfields_base_copy_comp      psc_mfields_cuda_copy_comp

#else
#error unknown FIELDS_BASE
#endif

#endif


