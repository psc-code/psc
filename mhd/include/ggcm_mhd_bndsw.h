
#ifndef GGCM_MHD_BND_SW_H
#define GGCM_MHD_BND_SW_H

#include <mrc_obj.h>

MRC_CLASS_DECLARE(ggcm_mhd_bndsw, struct ggcm_mhd_bndsw);

enum {
  SW_RR,
  SW_VX,
  SW_VY,
  SW_VZ,
  SW_PP,
  SW_BX,
  SW_BY,
  SW_BZ,
  SW_NR,
};

void ggcm_mhd_bndsw_new_step(struct ggcm_mhd_bndsw *bndsw);
void ggcm_mhd_bndsw_at(struct ggcm_mhd_bndsw *bndsw, float bntim, float xx[3],
		       float vals[SW_NR]);
void ggcm_mhd_bndsw_get_initial(struct ggcm_mhd_bndsw *bndsw,
				float vals[SW_NR]);

#endif
