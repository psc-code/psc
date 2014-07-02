
// #define BOUNDS_CHECK
#include <mrc_fld_as_double_aos.h>

#define CONDUCTING_Y2_OPS ggcm_mhd_bnd_conducting_y_ath_ops
#define CONDUCTING_Y2_FILL_GHOSTS ggcm_mhd_bnd_conducting_y_ath_fill_ghosts
#define CONDUCTING_Y2_STR "conducting_y_ath"
// #define CONDUCTING_Y2_SW 4

#include "ggcm_mhd_bnd_conducting_y2_common.c"
