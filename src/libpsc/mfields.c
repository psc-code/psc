
#include "psc.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

#define VAR(x) (void *)offsetof(struct psc_mfields_c, x)

static struct param psc_mfields_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  {},
};

#undef VAR

#define __MAKE_MFIELDS_METHODS(type)					\
									\
void									\
psc_mfields_##type##_set_domain(mfields_##type##_t *flds,		\
				struct mrc_domain *domain)		\
{									\
  psc_mfields_set_domain((struct psc_mfields *) flds, domain);		\
}									\
									\
void									\
psc_mfields_##type##_zero(mfields_##type##_t *flds, int m)		\
{									\
  psc_mfields_zero((struct psc_mfields *) flds, m);			\
}									\
									\
void									\
psc_mfields_##type##_set_comp(mfields_##type##_t *flds, int m, double alpha) \
{									\
  psc_mfields_set_comp((struct psc_mfields *) flds, m, alpha);		\
}									\
 									\
void									\
psc_mfields_##type##_scale(mfields_##type##_t *flds, double alpha)	\
{									\
  psc_mfields_scale((struct psc_mfields *) flds, alpha);		\
}									\
									\
void									\
psc_mfields_##type##_copy_comp(mfields_##type##_t *to, int mto,		\
			       mfields_##type##_t *from, int mfrom)	\
{									\
  psc_mfields_copy_comp((struct psc_mfields *) to, mto,			\
			(struct psc_mfields *) from, mfrom);		\
}									\
									\
void									\
psc_mfields_##type##_axpy(mfields_##type##_t *yf, double alpha,		\
			  mfields_##type##_t *xf)			\
{									\
  psc_mfields_axpy((struct psc_mfields *) yf, alpha,			\
		   (struct psc_mfields *) xf);				\
}									\
									\
LIST_HEAD(psc_mfields_##type##_list);					\
									\
void									\
psc_mfields_##type##_list_add(mfields_##type##_t **flds_p)		\
{									\
  psc_mfields_list_add(&psc_mfields_##type##_list,			\
		       (struct psc_mfields **) flds_p);			\
}									\
									\
void									\
psc_mfields_##type##_list_del(mfields_##type##_t **flds_p)		\
{									\
  psc_mfields_list_del(&psc_mfields_##type##_list,			\
		       (struct psc_mfields **) flds_p);			\
}									\


#define MAKE_MFIELDS_METHODS(type)					\
									\
static void								\
_psc_mfields_##type##_setup(mfields_##type##_t *flds)			\
{									\
  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,	\
						     &flds->nr_patches); \
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));			\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };	\
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],			\
		   patches[p].ldims[1] + flds->ibn[1],			\
		   patches[p].ldims[2] + flds->ibn[2] };		\
    fields_##type##_alloc(&flds->f[p], ilg, ihg, flds->nr_fields);	\
  }									\
}									\
									\
static void								\
_psc_mfields_##type##_destroy(mfields_##type##_t *flds)		        \
{									\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    fields_##type##_free(&flds->f[p]);					\
  }									\
  free(flds->f);							\
}									\
									\
struct mrc_class_psc_mfields_##type mrc_class_psc_mfields_##type = {	\
  .name             = "psc_mfields_" #type,				\
  .size             = sizeof(struct psc_mfields_##type),		\
  .param_descr      = psc_mfields_descr,				\
  .setup            = _psc_mfields_##type##_setup,			\
  .destroy          = _psc_mfields_##type##_destroy,			\
};

__MAKE_MFIELDS_METHODS(c)
__MAKE_MFIELDS_METHODS(fortran)
__MAKE_MFIELDS_METHODS(cuda)

MAKE_MFIELDS_METHODS(fortran)
//MAKE_MFIELDS_METHODS(sse2)

// ======================================================================

#define psc_mfields_ops(flds) (struct psc_mfields_ops *) ((flds)->obj.ops)

void
psc_mfields_set_domain(struct psc_mfields *flds, struct mrc_domain *domain)
{
  flds->domain = domain;
}

void
psc_mfields_zero(struct psc_mfields *flds, int m)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);
  assert(ops && ops->zero_comp);
  return ops->zero_comp(flds, m);
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
}

struct mrc_class_psc_mfields mrc_class_psc_mfields = {
  .name             = "psc_mfields",
  .size             = sizeof(struct psc_mfields),
  .init             = psc_mfields_init,
  .param_descr      = psc_mfields_descr,
  //  .setup            = _psc_mfields_setup,
  //  .destroy          = _psc_mfields_destroy,
};

