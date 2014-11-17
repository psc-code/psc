
#ifndef GGCM_MHD_CRDS_H
#define GGCM_MHD_CRDS_H

#include <mrc_obj.h>

#include <mrc_domain.h>

MRC_CLASS_DECLARE(ggcm_mhd_crds, struct ggcm_mhd_crds);

// cell-centered + other pre-calc coords as arrays
float *ggcm_mhd_crds_get_crd(struct ggcm_mhd_crds *crds, int d, int m);
float *ggcm_mhd_crds_get_global_crd(struct ggcm_mhd_crds *crds, int d);

// direct access to coords for a given cell
// (ideally avoided for performance critical parts, because it's slower)
void ggcm_mhd_crds_get_cc(struct ggcm_mhd_crds *crds, int ix, int iy, int iz, float crd_cc[3]);
void ggcm_mhd_crds_get_nc(struct ggcm_mhd_crds *crds, int ix, int iy, int iz, float crd_nc[3]);
void ggcm_mhd_crds_get_fc(struct ggcm_mhd_crds *crds, int ix, int iy, int iz, int d, float crd_fc[3]);
void ggcm_mhd_crds_get_ec(struct ggcm_mhd_crds *crds, int ix, int iy, int iz, int d, float crd_ec[3]);

#endif
