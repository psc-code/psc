
#ifndef GGCM_DIPOLE_MHD_H
#define GGCM_DIPOLE_MHD_H

#include <mrc_obj.h>

#include <mrc_fld.h>

// ======================================================================
// ggcm_mhd_dipole
//
// This object is responsible for setting up a dipolar magnetic fields
// in MHD fields

MRC_CLASS_DECLARE(ggcm_mhd_dipole, struct ggcm_mhd_dipole);

double ggcm_mhd_dipole_vector_potential(struct ggcm_mhd_dipole *mhd_dipole, int m,
					double x[3], float x0[3], float moment[3], float xmir);
void ggcm_mhd_dipole_add_dipole(struct ggcm_mhd_dipole *mhd_dipole, struct mrc_fld *b,
				float x0[3], float moment[3], float xmir, float keep);
void ggcm_mhd_dipole_set_b_field(struct ggcm_mhd_dipole *mhd_dipole,
				 float moment[3], double diptime);
void ggcm_mhd_dipole_update_b_field(struct ggcm_mhd_dipole *mhd_dipole,
				    struct mrc_fld *fld, double dacttime);

#endif
