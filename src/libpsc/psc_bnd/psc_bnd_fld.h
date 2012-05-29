
#include "psc.h"

void psc_bnd_fld_c_create(struct psc_bnd *bnd);
void psc_bnd_fld_c_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			      int mb, int me);
void psc_bnd_fld_c_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			       int mb, int me);

