
#ifndef PSC_SORT_H
#define PSC_SORT_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_sort, struct psc_sort);
void psc_sort_run(struct psc_sort *sort, struct psc_mparticles *particles);

END_C_DECLS

#endif
