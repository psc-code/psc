
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_collision;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_collision, struct psc_collision)

void psc_collision_run(struct psc_collision *collision, mparticles_base_t *particles);
