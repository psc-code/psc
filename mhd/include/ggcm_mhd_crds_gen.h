
#ifndef GGCM_MHD_CRDS_GEN_H
#define GGCM_MHD_CRDS_GEN_H

#include <mrc_obj.h>

#include "ggcm_mhd_crds.h"

MRC_CLASS_DECLARE(ggcm_mhd_crds_gen, struct ggcm_mhd_crds_gen);

void ggcm_mhd_crds_gen_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds);

#endif
