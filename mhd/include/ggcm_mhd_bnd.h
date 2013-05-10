
#ifndef GGCM_MHD_BND_H
#define GGCM_MHD_BND_H

#include <mrc_obj.h>

MRC_CLASS_DECLARE(ggcm_mhd_bnd, struct ggcm_mhd_bnd);

void ggcm_mhd_bnd_fill_ghosts(struct ggcm_mhd_bnd *bnd, int m, float bntim);

#endif
