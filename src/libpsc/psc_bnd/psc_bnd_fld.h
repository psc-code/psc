
#ifndef PSC_BND_FLD_H
#define PSC_BND_FLD_H

#include "psc.h"

BEGIN_C_DECLS

void psc_bnd_fld_c_create(struct psc_bnd *bnd);
void psc_bnd_fld_c_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
			      int mb, int me);
void psc_bnd_fld_c_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
			       int mb, int me);

void psc_bnd_fld_single_create(struct psc_bnd *bnd);
void psc_bnd_fld_single_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
				   int mb, int me);
void psc_bnd_fld_single_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
				    int mb, int me);

void psc_bnd_fld_mix_create(struct psc_bnd *bnd);
void psc_bnd_fld_mix_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
				int mb, int me);
void psc_bnd_fld_mix_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base,
				 int mb, int me);


void psc_bnd_fld_single_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);
void psc_bnd_fld_single_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);
void psc_bnd_fld_single_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);

void psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);
void psc_bnd_fld_cuda_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);
void psc_bnd_fld_cuda_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx);

END_C_DECLS

#endif
