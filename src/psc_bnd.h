
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_bnd;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_bnd, struct psc_bnd)

void psc_bnd_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me);
void psc_bnd_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me);
void psc_bnd_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles);
