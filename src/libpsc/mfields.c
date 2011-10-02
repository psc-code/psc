
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
  flds->domain = domain;						\
}									\
									\
static inline struct psc_mfields_##type##_ops *				\
to_##type##_ops(mfields_##type##_t *flds)				\
{									\
  return (struct psc_mfields_##type##_ops *) flds->obj.ops;		\
}									\
									\
void									\
psc_mfields_##type##_zero(mfields_##type##_t *flds, int m)		\
{									\
  struct psc_mfields_##type##_ops *ops = to_##type##_ops(flds);		\
  assert(ops && ops->zero_comp);					\
  return ops->zero_comp(flds, m);					\
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
  .init             = psc_mfields_##type##_init,			\
  .param_descr      = psc_mfields_descr,				\
  .setup            = _psc_mfields_##type##_setup,			\
  .destroy          = _psc_mfields_##type##_destroy,			\
};

__MAKE_MFIELDS_METHODS(c)

static void
psc_mfields_fortran_init()
{
}

__MAKE_MFIELDS_METHODS(fortran)
__MAKE_MFIELDS_METHODS(cuda)

MAKE_MFIELDS_METHODS(fortran)
//MAKE_MFIELDS_METHODS(sse2)

