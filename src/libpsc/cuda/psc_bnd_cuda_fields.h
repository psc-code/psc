
#ifndef PSC_BND_CUDA_FIELDS_H
#define PSC_BND_CUDA_FIELDS_H

#include "psc_bnd.h"

void psc_bnd_fld_cuda_create(struct psc_bnd *bnd);
void psc_bnd_fld_cuda_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me);
void psc_bnd_fld_cuda_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me);
void psc_bnd_fld_cuda_add_ghosts_prep(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
void psc_bnd_fld_cuda_add_ghosts_post(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
void psc_bnd_fld_cuda_fill_ghosts_prep(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
void psc_bnd_fld_cuda_fill_ghosts_post(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);

#endif

