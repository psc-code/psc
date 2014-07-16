
#ifndef MHD_RIEMANN_H
#define MHD_RIEMANN_H

#include <mrc_fld.h>

// ======================================================================
// mhd_riemann

MRC_CLASS_DECLARE(mhd_riemann, struct mhd_riemann);

void
mhd_riemann_run(struct mhd_riemann *riemann,
		struct mrc_fld *F1d,
		struct mrc_fld *Ul, struct mrc_fld *Ur,
		struct mrc_fld *Wl, struct mrc_fld *Wr,
		int ldim, int l, int r, int dir);

#endif
