
#include <mrc_fld_as_double_aos.h>
#define MT MT_FULLY_CONSERVATIVE

#define ggcm_mhd_bnd_ops_sphere ggcm_mhd_bnd_ops_sphere_fc_double_aos
#define ggcm_mhd_bnd_sub_name "sphere_fc_double_aos"

#include "ggcm_mhd_bnd_sphere_common.c"
