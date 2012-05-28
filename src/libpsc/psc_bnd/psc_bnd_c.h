
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_ddc.h>

void psc_bnd_fields_c_create(struct psc_bnd *bnd);
void psc_bnd_fields_c_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				 int mb, int me);
void psc_bnd_fields_c_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				  int mb, int me);

