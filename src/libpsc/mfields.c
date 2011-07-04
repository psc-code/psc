
#include "psc.h"

#include <stdlib.h>

#define MAKE_MFIELDS_METHODS(type)					\
									\
LIST_HEAD(mfields_##type##_list);					\
									\
void									\
psc_mfields_##type##_set_domain(mfields_##type##_t *flds, struct mrc_domain *domain, \
			    int nr_fields, int ibn[3])			\
{									\
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &flds->nr_patches); \
  flds->f = calloc(flds->nr_patches, sizeof(*flds->f));			\
  for (int p = 0; p < flds->nr_patches; p++) {				\
    int ilg[3] = { -ibn[0], -ibn[1], -ibn[2] };				\
    int ihg[3] = { patches[p].ldims[0] + ibn[0],			\
		   patches[p].ldims[1] + ibn[1],			\
		   patches[p].ldims[2] + ibn[2] };			\
    fields_##type##_alloc(&flds->f[p], ilg, ihg, nr_fields);		\
  }									\
  list_add_tail(&flds->entry, &mfields_##type##_list);			\
}									\
									\
static void									\
_psc_mfields_##type##_destroy(mfields_##type##_t *flds)				\
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
  .destroy          = _psc_mfields_##type##_destroy,			\
};

MAKE_MFIELDS_METHODS(fortran)
MAKE_MFIELDS_METHODS(c)
//MAKE_MFIELDS_METHODS(sse2)
