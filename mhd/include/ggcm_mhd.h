
#ifndef GGCM_MHD_H
#define GGCM_MHD_H

#include <mrc_obj.h>

#include "ggcm_mhd_flds.h"

// ======================================================================
// ggcm_mhd
//
// This object runs an MHD simulation

MRC_CLASS_DECLARE(ggcm_mhd, struct ggcm_mhd);

void ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, int m, float bntim);
void ggcm_mhd_newstep(struct ggcm_mhd *mhd, float *dtn);
void ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct ggcm_mhd_flds *flds,
			struct mrc_f3 *divg);
void ggcm_mhd_calc_currcc(struct ggcm_mhd *mhd, struct ggcm_mhd_flds *flds, int m,
			struct mrc_f3 *currcc);
void ggcm_mhd_push(struct ggcm_mhd *mhd, float *dtn,
		   bool do_nwst, bool do_iono, bool do_rcm);
void ggcm_mhd_push_step(struct ggcm_mhd *mhd);
void ggcm_mhd_get_state(struct ggcm_mhd *mhd);
void ggcm_mhd_set_state(struct ggcm_mhd *mhd);

int ggcm_mhd_ntot(struct ggcm_mhd *mhd);

void ts_ggcm_mhd_step_calc_rhs(void *ctx, struct mrc_obj *_rhs, float time,
			       struct mrc_obj *_fld);

// ----------------------------------------------------------------------

void ggcm_mhd_register();

#endif
