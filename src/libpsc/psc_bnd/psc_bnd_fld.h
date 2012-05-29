
#ifndef PSC_BND_FLD_H
#define PSC_BND_FLD_H

#include "psc.h"

void psc_bnd_fld_c_create(struct psc_bnd *bnd);
void psc_bnd_fld_c_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			      int mb, int me);
void psc_bnd_fld_c_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
			       int mb, int me);

void psc_bnd_fld_single_create(struct psc_bnd *bnd);
void psc_bnd_fld_single_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				   int mb, int me);
void psc_bnd_fld_single_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				    int mb, int me);

void psc_bnd_fld_mix_create(struct psc_bnd *bnd);
void psc_bnd_fld_mix_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				int mb, int me);
void psc_bnd_fld_mix_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base,
				 int mb, int me);

#endif
