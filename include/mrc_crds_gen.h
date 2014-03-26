
#ifndef MRC_CRDS_GEN_H
#define MRC_CRDS_GEN_H

#include <mrc_obj.h>

MRC_CLASS_DECLARE(mrc_crds_gen, struct mrc_crds_gen);

void mrc_crds_gen_run(struct mrc_crds_gen *gen, double *xx, double *dx);

#endif

