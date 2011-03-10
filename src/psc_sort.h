
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_sort;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_sort, struct psc_sort)

void psc_sort_run(struct psc_sort *sort, mparticles_base_t *particles);
