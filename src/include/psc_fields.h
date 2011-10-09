
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
struct psc_mfields *_psc_mfields_c_get_fortran(struct psc_mfields *_flds_base, int mb, int me);
void _psc_mfields_c_put_fortran(struct psc_mfields *flds, struct psc_mfields *_flds_base, int mb, int me);

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

extern struct psc_mfields_ops psc_mfields_c_ops;
extern struct psc_mfields_ops psc_mfields_fortran_ops;
extern struct psc_mfields_ops psc_mfields_cuda_ops;

typedef struct psc_mfields mfields_fortran_t;
static inline fields_fortran_t *
psc_mfields_get_patch_fortran(struct psc_mfields *flds, int p)
{
  assert((void *) psc_mfields_ops(flds) == (void *) &psc_mfields_fortran_ops);
  return ((fields_fortran_t *)flds->data) + p;
}

#include "psc_fields_c.h"
typedef struct psc_mfields mfields_c_t;
static inline fields_c_t *
psc_mfields_get_patch_c(struct psc_mfields *flds, int p)
{
  assert((void *) psc_mfields_ops(flds) == (void *) &psc_mfields_c_ops);
  return ((fields_c_t *)flds->data) + p;
}

#ifdef USE_SSE2
#include "psc_fields_sse2.h"
typedef struct psc_mfields mfields_sse2_t;
#endif

#ifdef USE_CUDA
#include "psc_fields_cuda.h"
typedef struct psc_mfields mfields_cuda_t;
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


