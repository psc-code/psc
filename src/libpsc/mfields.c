
#include "psc.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================
// _psc_mfields_setup

static void
_psc_mfields_setup(struct psc_mfields *flds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);

  flds->comp_name = calloc(flds->nr_fields, sizeof(*flds->comp_name));

  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->flds = calloc(flds->nr_patches, sizeof(*flds->flds));
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = psc_fields_create(psc_mfields_comm(flds));
    psc_fields_set_type(pf, ops->name);
    for (int d = 0; d < 3; d++) {
      pf->ib[d] = -flds->ibn[d];
      pf->im[d] = patches[p].ldims[d] + 2 * flds->ibn[d];
    }
    pf->nr_comp = flds->nr_fields;
    pf->first_comp = flds->first_comp;
    psc_fields_setup(pf);
    flds->flds[p] = pf;
  }
}

// ======================================================================
// psc_mfields_destroy

static void
_psc_mfields_destroy(struct psc_mfields *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = flds->flds[p];
    psc_fields_destroy(pf);
  }
  free(flds->flds);

  // sub-destroy has already been called
  for (int m = 0; m < flds->nr_fields; m++) {
    free(flds->comp_name[m]);
  }
  free(flds->comp_name);
}

void
psc_mfields_set_domain(struct psc_mfields *flds, struct mrc_domain *domain)
{
  flds->domain = domain;
}

void
psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s)
{
  assert(m >= flds->first_comp && m < flds->first_comp + flds->nr_fields);
  flds->comp_name[m - flds->first_comp] = strdup(s);
}

const char *
psc_mfields_comp_name(struct psc_mfields *flds, int m)
{
  assert(m >= flds->first_comp && m < flds->first_comp + flds->nr_fields);
  return flds->comp_name[m - flds->first_comp];
}

void
psc_mfields_zero(struct psc_mfields *flds, int m)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);
  assert(ops && ops->zero_comp);
  return ops->zero_comp(flds, m);
}

void
psc_mfields_zero_range(struct psc_mfields *flds, int mb, int me)
{
  for (int m = mb; m < me; m++) {
    psc_mfields_zero(flds, m);
  }
}

void
psc_mfields_set_comp(struct psc_mfields *flds, int m, double alpha)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);
  assert(ops && ops->set_comp);
  return ops->set_comp(flds, m, alpha);
}

void
psc_mfields_scale(struct psc_mfields *flds, double alpha)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);
  assert(ops && ops->scale);
  return ops->scale(flds, alpha);
}

void
psc_mfields_copy_comp(struct psc_mfields *to, int mto,
		      struct psc_mfields *from, int mfrom)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(to);
  assert(ops == psc_mfields_ops(from));
  assert(ops && ops->copy_comp);
  return ops->copy_comp(to, mto, from, mfrom);
}

void
psc_mfields_axpy(struct psc_mfields *yf, double alpha,
		 struct psc_mfields *xf)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(yf);
  assert(ops == psc_mfields_ops(xf));
  assert(ops && ops->axpy);
  return ops->axpy(yf, alpha, xf);
}

void
psc_mfields_axpy_comp(struct psc_mfields *yf, int ym, double alpha,
		      struct psc_mfields *xf, int xm)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(yf);
  assert(ops == psc_mfields_ops(xf));
  assert(ops && ops->axpy_comp);
  return ops->axpy_comp(yf, ym, alpha, xf, xm);
}

#define MAKE_MFIELDS_GET_PUT(type)					\
									\
struct psc_mfields *							\
psc_mfields_get_##type(struct psc_mfields *base, int mb, int me)	\
{									\
  struct psc_mfields_ops *ops = psc_mfields_ops(base);			\
  if (ops == &psc_mfields_##type##_ops) {				\
    return base;							\
  }									\
  assert(ops);								\
  static int pr;							\
  if (!pr) {								\
    pr = prof_register("mfields_get_" #type, 1., 0, 0);			\
  }									\
  prof_start(pr);							\
  struct psc_mfields *flds = psc_mfields_create(psc_mfields_comm(base)); \
  psc_mfields_set_type(flds, #type);					\
  psc_mfields_set_domain(flds, base->domain);				\
  psc_mfields_set_param_int(flds, "nr_fields", base->nr_fields);	\
  psc_mfields_set_param_int3(flds, "ibn", base->ibn);			\
  psc_mfields_set_param_int(flds, "first_comp", base->first_comp);	\
  psc_mfields_setup(flds);						\
									\
  if (!ops->copy_to_##type) {						\
    fprintf(stderr, "ERROR: missing copy_to_"#type" in psc_mfields '%s'\n", \
	    psc_mfields_type(base));					\
    assert(0);								\
  }									\
  ops->copy_to_##type(base, flds, mb, me);				\
  prof_stop(pr);							\
									\
  return flds;								\
}									\
									\
void									\
psc_mfields_put_##type(struct psc_mfields *flds,			\
		       struct psc_mfields *base, int mb, int me)	\
{									\
  struct psc_mfields_ops *ops = psc_mfields_ops(base);			\
  struct psc_mfields_ops *ops2 = psc_mfields_ops(flds);			\
  if (ops == ops2) {							\
    return;								\
  }									\
  assert(ops && ops2);							\
  static int pr;							\
  if (!pr) {								\
    pr = prof_register("mfields_put_" #type, 1., 0, 0);			\
  }									\
  prof_start(pr);							\
									\
  if (!ops->copy_from_##type) {						\
    fprintf(stderr, "ERROR: missing copy_from_"#type" in psc_mfields '%s'!\n", \
	    psc_mfields_type(base));					\
    assert(0);								\
  }									\
  ops->copy_from_##type(base, flds, mb, me);				\
  psc_mfields_destroy(flds);						\
  prof_stop(pr);							\
}									\

MAKE_MFIELDS_GET_PUT(c)
MAKE_MFIELDS_GET_PUT(fortran)
#ifdef USE_CUDA
MAKE_MFIELDS_GET_PUT(cuda)
#endif

// ======================================================================

void
psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p)
{
  struct psc_mfields_list_entry *p = malloc(sizeof(*p));
  p->flds_p = flds_p;
  list_add_tail(&p->entry, head);
}

void
psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p)
{
  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, head, entry, struct psc_mfields_list_entry) {
    if (p->flds_p == flds_p) {
      list_del(&p->entry);
      free(p);
      return;
    }
  }
  assert(0);
}

// ======================================================================

static void
psc_mfields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_cuda_ops);
#endif
}

#define VAR(x) (void *)offsetof(struct psc_mfields, x)
static struct param psc_mfields_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  { "first_comp"     , VAR(first_comp)      , PARAM_INT(0)        },
  {},
};
#undef VAR

struct mrc_class_psc_mfields mrc_class_psc_mfields = {
  .name             = "psc_mfields",
  .size             = sizeof(struct psc_mfields),
  .init             = psc_mfields_init,
  .param_descr      = psc_mfields_descr,
  .setup            = _psc_mfields_setup,
  .destroy          = _psc_mfields_destroy,
};

