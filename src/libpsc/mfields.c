
#include "psc.h"

#include <mrc_params.h>
#include <stdlib.h>

#define VAR(x) (void *)offsetof(struct psc_mfields_c, x)

static struct param psc_mfields_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  {},
};

#undef VAR

#define MAKE_MFIELDS_METHODS(type)					\
									\
LIST_HEAD(mfields_##type##_list);					\
									\
void									\
psc_mfields_##type##_set_domain(mfields_##type##_t *flds,		\
				struct mrc_domain *domain)		\
{									\
  flds->domain = domain;						\
}									\
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
  list_add_tail(&flds->entry, &mfields_##type##_list);			\
}									\
									\
static void								\
_psc_mfields_##type##_destroy(mfields_##type##_t *flds)		\
{									\
  if (!flds)								\
    return;								\
									\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    fields_##type##_free(&flds->f[p]);					\
  }									\
  free(flds->f);							\
  list_del(&flds->entry);						\
  free(flds);								\
}									\
									\
struct mrc_class_psc_mfields_##type mrc_class_psc_mfields_##type = {	\
  .name             = "psc_mfields_" #type,				\
  .size             = sizeof(struct psc_mfields_##type),		\
  .param_descr      = psc_mfields_descr,				\
  .setup            = _psc_mfields_##type##_setup,			\
  .destroy          = _psc_mfields_##type##_destroy,			\
};

MAKE_MFIELDS_METHODS(fortran)
MAKE_MFIELDS_METHODS(c)
//MAKE_MFIELDS_METHODS(sse2)
