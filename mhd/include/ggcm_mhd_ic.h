
#ifndef GGCM_MHD_IC_H
#define GGCM_MHD_IC_H

#include <mrc_obj.h>

#include "ggcm_mhd.h"

MRC_CLASS_DECLARE(ggcm_mhd_ic, struct ggcm_mhd_ic);

void ggcm_mhd_ic_set_mhd(struct ggcm_mhd_ic *ic, struct ggcm_mhd *mhd);
void ggcm_mhd_ic_run(struct ggcm_mhd_ic *ic);

#endif
