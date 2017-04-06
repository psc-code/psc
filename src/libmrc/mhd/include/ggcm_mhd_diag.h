
#ifndef GGCM_MHD_DIAG_H
#define GGCM_MHD_DIAG_H

#include <mrc_obj.h>

#include <mrc_fld.h>
#include <mrc_mod.h>

// ======================================================================
// ggcm_mhd_diag
//
// This object is responsible for handling output from the MHD nodes

MRC_CLASS_DECLARE(ggcm_mhd_diag, struct ggcm_mhd_diag);

// Output whatever is needed at this time.
// (called at every timestep, most of the time it doesn't
// actually write any output, but just returns.)
void ggcm_mhd_diag_run(struct ggcm_mhd_diag *diag);

// Run output of the specified field (which should be an MHD state) now,
// of the specific type diag_type and with an output index of itdia.
void ggcm_mhd_diag_run_now(struct ggcm_mhd_diag *diag, struct mrc_fld *fld,
			   int diag_type, int itdia);

// Shutdown the output server (if any).
void ggcm_mhd_diag_shutdown(struct ggcm_mhd_diag *diag);

// register mrc_mod for actual diag server
void ggcm_mhd_diag_mod_register(struct ggcm_mhd_diag *mhd_diag, struct mrc_mod *mod);

#endif
