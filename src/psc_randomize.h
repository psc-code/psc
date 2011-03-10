
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_randomize;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_randomize, struct psc_randomize)

void psc_randomize_run(struct psc_randomize *randomize, mparticles_base_t *particles);
