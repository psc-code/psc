
#ifndef GGCM_MHD_STEP_H
#define GGCM_MHD_STEP_H

#include <mrc_obj.h>

#include <mrc_fld.h>

// ======================================================================
// ggcm_mhd_step
//
// This class is responsible for doing a single MHD timestep

MRC_CLASS_DECLARE(ggcm_mhd_step, struct ggcm_mhd_step);

// calculate the r.h.s. in the given scheme
// should be more general than taking mrc_fld, eventually
void ggcm_mhd_step_calc_rhs(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
			    struct mrc_fld *x);

// perform one entire time step in the given scheme
void ggcm_mhd_step_run(struct ggcm_mhd_step *step, struct mrc_fld *x);

int ggcm_mhd_step_mhd_type(struct ggcm_mhd_step *step);

#endif
