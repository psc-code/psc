
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_ddc.h>

struct mrc_ddc *psc_bnd_lib_create_ddc(struct psc *psc);
void __psc_bnd_lib_add_ghosts(struct mrc_ddc *ddc, mfields_t *flds,
			      int mb, int me);
void psc_bnd_lib_add_ghosts(struct mrc_ddc *ddc, mfields_base_t *flds_base,
			    int mb, int me);
void psc_bnd_lib_fill_ghosts(struct mrc_ddc *ddc, mfields_base_t *flds_base,
			     int mb, int me);

