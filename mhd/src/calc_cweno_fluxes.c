#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <assert.h>
#include <stdio.h>
#include <math.h>

  
// ----------------------------------------------------------------------
// calc_cweno_fluxes
//
// calculates CWENO fluxes on faces in flux_E, from the original state
// vector u (which is cell centered / on the Yee grid)
// flux[0-4] are the fluid vars
// flux[5-7] are E-field "fluxes" (not B-field!)
// (Ziegler 2004 section 3.1) 

void
calc_cweno_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
		  struct mrc_fld *_u)
{
  //initialize cell surface center variables			    
  struct mrc_fld *_u_delta[3], *_u_p[3], *_u_m[3], *_E_p[3], *_E_m[3];
  struct mrc_fld *_flux_p[3], *_flux_m[3];
  for (int f = 0; f < 3; f++) {
    _u_delta[f] = ggcm_mhd_get_fields(mhd, "u_delta", _B1Z + 1);
    _u_p[f] = ggcm_mhd_get_fields(mhd, "u_p", _JZ + 1);
    _u_m[f] = ggcm_mhd_get_fields(mhd, "u_m", _JZ + 1);
    _E_p[f] = ggcm_mhd_get_fields(mhd, "E_p", 3);
    _E_m[f] = ggcm_mhd_get_fields(mhd, "E_m", 3);
    _flux_p[f] = ggcm_mhd_get_fields(mhd, "flux_p", 8);
    _flux_m[f] = ggcm_mhd_get_fields(mhd, "flux_m", 8);
  }

  ggcm_mhd_fill_ghosts(mhd, _u, 0, mhd->time);

  calc_u_delta(_u_delta, _u);
  calc_u_pm(mhd, _u_p, _u_m, _E_p, _E_m, _u, _u_delta);
  
#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "u_p0");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, _u_p[0], DIAG_TYPE_3D, cnt++);
  }
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "u_m0");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, _u_m[0], DIAG_TYPE_3D, cnt++);
  }
#endif

  // calculate fluxes per face (small f's) using reconstructed 
  // variables U^N(SWETB) and B^N(SWETB) = (Bx,By,Bz)^N(SEWTB)
  for (int f = 0; f < 3; f++) {
    calc_fluxes_per_face(_flux_p, mhd, _u_p[f], f);
    calc_fluxes_per_face(_flux_m, mhd, _u_m[f], f);
  }
  
  calc_KNP_fluxes(mhd, _flux, _flux_p, _flux_m, _u, _u_p, _u_m, _E_p, _E_m);

  for (int f = 0; f < 3; f++) {
    mrc_fld_destroy(_flux_p[f]);
    mrc_fld_destroy(_flux_m[f]);
    mrc_fld_destroy(_u_delta[f]);
    mrc_fld_destroy(_u_p[f]);
    mrc_fld_destroy(_u_m[f]);
    mrc_fld_destroy(_E_p[f]);
    mrc_fld_destroy(_E_m[f]);
  }
}
