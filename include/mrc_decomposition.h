
#ifndef MRC_DECOMPOSITION_H
#define MRC_DECOMPOSITION_H

#include <mrc_obj.h>

// ======================================================================
// mrc_decomposition

MRC_CLASS_DECLARE(mrc_decomposition, struct mrc_decomposition);

int mrc_decomposition_global_to_local(struct mrc_decomposition *dc, int gidx);
bool mrc_decomposition_is_local(struct mrc_decomposition *dc, int gidx);
int mrc_decomposition_find_rank(struct mrc_decomposition *dc, int gidx);

#endif
