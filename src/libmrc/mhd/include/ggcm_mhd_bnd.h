
#ifndef GGCM_MHD_BND_H
#define GGCM_MHD_BND_H

#include <mrc_obj.h>
#include <mrc_fld.h>

MRC_CLASS_DECLARE(ggcm_mhd_bnd, struct ggcm_mhd_bnd);

void ggcm_mhd_bnd_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
			      int m, float bntim);
void ggcm_mhd_bnd_fill_ghosts_E(struct ggcm_mhd_bnd *bnd, struct mrc_fld *E);

#endif
