
#ifndef GGCM_MHD_STEP_H
#define GGCM_MHD_STEP_H

#include <mrc_obj.h>

#include <mrc_fld.h>

// ======================================================================
// ggcm_mhd_step
//
// This class is responsible for doing a single MHD timestep

MRC_CLASS_DECLARE(ggcm_mhd_step, struct ggcm_mhd_step);

double ggcm_mhd_step_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x);

// calculate the r.h.s. in the given scheme
// should be more general than taking mrc_fld, eventually
void ggcm_mhd_step_calc_rhs(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
			    struct mrc_fld *x);

// calculate electric field (edge centered) for daig, since each scheme
// treats E differently. This will probably have to set things up, so it's
// best not to call this every time step.
void ggcm_mhd_step_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *E,
                            struct mrc_fld *x);

// perform one entire time step in the given scheme
void ggcm_mhd_step_run(struct ggcm_mhd_step *step, struct mrc_fld *x);

// sets up mhd->fld and aux fields as needed by the selected step subclass
void ggcm_mhd_step_setup_flds(struct ggcm_mhd_step *step);

// returns whether this particular ggcm_mhd_step implementation provides just a
// calc_rhs() function (needs a timestepper like "rk2"), or provides a run() function
// that takes over the entire step (then the timestepper should be just "step").
bool ggcm_mhd_step_has_calc_rhs(struct ggcm_mhd_step *step);

#endif
