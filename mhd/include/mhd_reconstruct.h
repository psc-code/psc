
#ifndef MHD_RECONSTRUCT_H
#define MHD_RECONSTRUCT_H

#include <mrc_obj.h>

#include <mrc_fld.h>

// ======================================================================
// mhd_reconstruct

MRC_CLASS_DECLARE(mhd_reconstruct, struct mhd_reconstruct);

void mhd_reconstruct_run(struct mhd_reconstruct *mr,
			 struct mrc_fld *Ul, struct mrc_fld *Ur,
			 struct mrc_fld *Wl, struct mrc_fld *Wr,
			 struct mrc_fld *W1d, struct mrc_fld *Bxi,
			 int ldim, int l, int r, int dir);

#endif
