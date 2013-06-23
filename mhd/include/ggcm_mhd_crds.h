
#ifndef GGCM_MHD_CRDS_H
#define GGCM_MHD_CRDS_H

#include <mrc_obj.h>

#include <mrc_domain.h>

MRC_CLASS_DECLARE(ggcm_mhd_crds, struct ggcm_mhd_crds);

float *ggcm_mhd_crds_get_crd(struct ggcm_mhd_crds *crds, int d, int m);

#endif
