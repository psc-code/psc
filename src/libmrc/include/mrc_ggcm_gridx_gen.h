
#ifndef MRC_GGCM_GRIDX_GEN_H
#define MRC_GGCM_GRIDX_GEN_H

#include "mrc_crds_gen_private.h"

void
generate_ggcm_x_grid(struct mrc_crds_gen *gen, double *xx, double *dx,
                     double (*dx_func)(struct mrc_crds_gen *gen, double x, double fak));

#endif
