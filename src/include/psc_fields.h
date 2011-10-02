
#ifndef PSC_FIELDS_H
#define PSC_FIELDS_H

// ----------------------------------------------------------------------
// mfields type

struct psc_mfields {
  struct mrc_obj obj;
  void *data;
  int nr_patches;
  struct mrc_domain *domain;
  int nr_fields;
  int ibn[3];
};

MRC_CLASS_DECLARE(psc_mfields, struct psc_mfields);

struct psc_mfields_ops {
  MRC_SUBCLASS_OPS(struct psc_mfields);
  void (*zero_comp)(struct psc_mfields *, int m);
  void (*set_comp)(struct psc_mfields *, int m, double alpha);
  void (*scale)(struct psc_mfields *, double alpha);
  void (*copy_comp)(struct psc_mfields *to, int mto,
		    struct psc_mfields *from, int mfrom);
  void (*axpy)(struct psc_mfields *yf, double alpha,
	       struct psc_mfields *xf);
  struct psc_mfields *(*get_c)(struct psc_mfields *base, int mb, int me);
  struct psc_mfields *(*get_fortran)(struct psc_mfields *base, int mb, int me);
  struct psc_mfields *(*get_cuda)(struct psc_mfields *base, int mb, int me);
  void (*put_c)(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me);
  void (*put_fortran)(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me);
  void (*put_cuda)(struct psc_mfields *flds, struct psc_mfields *flds_base, int mb, int me);
};

void psc_mfields_set_domain(struct psc_mfields *flds,
			    struct mrc_domain *domain);
void psc_mfields_zero(struct psc_mfields *flds, int m);
void psc_mfields_set_comp(struct psc_mfields *flds, int m, double alpha);
void psc_mfields_scale(struct psc_mfields *flds, double alpha);
void psc_mfields_copy_comp(struct psc_mfields *to, int mto,
			   struct psc_mfields *from, int mfrom);
void psc_mfields_axpy(struct psc_mfields *yf, double alpha,
		      struct psc_mfields *xf);
struct psc_mfields *psc_mfields_get_c(struct psc_mfields *base, int mb, int me);
struct psc_mfields *psc_mfields_get_fortran(struct psc_mfields *base, int mb, int me);
struct psc_mfields *psc_mfields_get_cuda(struct psc_mfields *base, int mb, int me);
void psc_mfields_put_c(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me);
void psc_mfields_put_fortran(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me);
void psc_mfields_put_cuda(struct psc_mfields *flds, struct psc_mfields *base, int mb, int me);

struct psc_mfields_list_entry {
  struct psc_mfields **flds_p;
  list_t entry;
};

void psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p);
void psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p);

struct psc_mfields *_psc_mfields_c_get_cuda(struct psc_mfields *_flds_base, int mb, int me);
void _psc_mfields_c_put_cuda(struct psc_mfields *flds, struct psc_mfields *_flds_base, int mb, int me);

// This type is replicated for each actual fields type, however,
// the interface and implementation is always identical, hence 
// created automatically for the variants using macros

#define DECLARE_MFIELDS_METHODS(type)					\
  struct psc_mfields_##type {						\
    struct mrc_obj obj;							\
    fields_##type##_t *xf;						\
    int nr_patches;							\
    struct mrc_domain *domain;						\
    int nr_fields;							\
    int ibn[3];								\
  };									\
  typedef struct psc_mfields mfields_##type##_t;			\
  									\
  MRC_CLASS_DECLARE(psc_mfields_##type, struct psc_mfields_##type);	\
									\
  struct psc_mfields_##type##_ops {					\
    MRC_SUBCLASS_OPS(struct psc_mfields_##type);			\
    void (*zero_comp)(struct psc_mfields_##type *, int m);		\
    void (*set_comp)(struct psc_mfields_##type *, int m, double alpha);	\
    void (*scale)(struct psc_mfields_##type *, double alpha);		\
    void (*copy_comp)(struct psc_mfields_##type *to, int mto,		\
		      struct psc_mfields_##type *from, int mfrom);	\
    void (*axpy)(struct psc_mfields_##type *yf, double alpha,		\
		 struct psc_mfields_##type *xf);			\
    struct psc_mfields *(*get_c)(struct psc_mfields *flds_base, int mb, int me); \
    struct psc_mfields *(*get_fortran)(struct psc_mfields *flds_base, int mb, int me); \
    struct psc_mfields *(*get_cuda)(struct psc_mfields *flds_base, int mb, int me); \
    void (*put_c)(struct psc_mfields *flds, struct psc_mfields *flds_base, int mb, int me); \
    void (*put_fortran)(struct psc_mfields *flds, struct psc_mfields *flds_base, int mb, int me); \
    void (*put_cuda)(struct psc_mfields *flds, struct psc_mfields *flds_base, int mb, int me); \
  };									\
  									\
  mfields_##type##_t *							\
  psc_mfields_##type##_get_from(int mb, int me, void *flds_base);	\
  void psc_mfields_##type##_put_to(mfields_##type##_t *pf, int mb, int me, void *flds_base); \

  

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

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

extern struct psc_mfields_ops psc_mfields_c_ops;
extern struct psc_mfields_ops psc_mfields_fortran_ops;
extern struct psc_mfields_ops psc_mfields_cuda_ops;

static inline fields_c_t *
psc_mfields_get_patch_c(struct psc_mfields *flds, int p)
{
  assert((void *) psc_mfields_ops(flds) == (void *) &psc_mfields_c_ops);
  return ((fields_c_t *)flds->data) + p;
}

static inline fields_fortran_t *
psc_mfields_get_patch_fortran(struct psc_mfields *flds, int p)
{
  assert((void *) psc_mfields_ops(flds) == (void *) &psc_mfields_fortran_ops);
  return ((fields_fortran_t *)flds->data) + p;
}

#ifdef USE_CUDA
static inline fields_cuda_t *
psc_mfields_get_patch_cuda(struct psc_mfields *flds, int p)
{
  assert((void *) psc_mfields_ops(flds) == (void *) &psc_mfields_cuda_ops);
  return ((fields_cuda_t *)flds->data) + p;
}

#endif

extern list_t psc_mfields_base_list;

// ----------------------------------------------------------------------
// base fields type

typedef struct psc_mfields mfields_base_t;

#if FIELDS_BASE == FIELDS_FORTRAN

#define s_fields_base "fortran"

#elif FIELDS_BASE == FIELDS_C

#define s_fields_base "c"

#elif FIELDS_BASE == FIELDS_SSE2

#include "psc_fields_sse2.h"

#elif FIELDS_BASE == FIELDS_CUDA

#include "psc_fields_cuda.h"
#define s_fields_base "cuda"

#else
#error unknown FIELDS_BASE
#endif

#endif


