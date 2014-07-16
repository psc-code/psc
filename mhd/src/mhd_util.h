
#ifndef MHD_UTIL_H
#define MHD_UTIL_H

#include "ggcm_mhd.h"

void setup_mrc_fld_1d(struct mrc_fld *f, struct mrc_fld *f_tmpl, int nr_comps);
void setup_mrc_fld_3d(struct mrc_fld *f, struct mrc_fld *f_tmpl, int nr_comps);

#endif
