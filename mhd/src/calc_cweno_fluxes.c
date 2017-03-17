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
calc_cweno_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *flux[3],
		  struct mrc_fld *u)
{
  //initialize cell surface center variables			    
  struct mrc_fld *u_delta[3], *u_p[3], *u_m[3], *E_p[3], *E_m[3];
  struct mrc_fld *flux_p[3], *flux_m[3];
  
  for (int f = 0; f < 3; f++) {
    u_delta[f] = ggcm_mhd_get_fields(mhd, "u_delta", BZ + 1);
    u_p[f] = ggcm_mhd_get_fields(mhd, "u_p", _JZ + 1);
    u_m[f] = ggcm_mhd_get_fields(mhd, "u_m", _JZ + 1);
    E_p[f] = ggcm_mhd_get_fields(mhd, "E_p", 3);
    E_m[f] = ggcm_mhd_get_fields(mhd, "E_m", 3);
    flux_p[f] = ggcm_mhd_get_fields(mhd, "flux_p", 8);
    flux_m[f] = ggcm_mhd_get_fields(mhd, "flux_m", 8);
  }
  
  ggcm_mhd_fill_ghosts(mhd, u, mhd->time_code);
  calc_u_delta(u_delta, u, mhd);
  
  ggcm_mhd_fill_ghosts(mhd, u_delta[0], mhd->time_code); 
  ggcm_mhd_fill_ghosts(mhd, u_delta[1], mhd->time_code); 
  ggcm_mhd_fill_ghosts(mhd, u_delta[2], mhd->time_code); 
  
#if CWENOREC
  calc_u_cweno(mhd, u_p, u_m, E_p, E_m, u, u_delta);
#else 
  calc_u_pm(mhd, u_p, u_m, E_p, E_m, u, u_delta);
#endif 
  
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
    calc_fluxes_per_face(flux_p, mhd, u_p[f], f);
    calc_fluxes_per_face(flux_m, mhd, u_m[f], f);
  }
  
  calc_KNP_fluxes(mhd, flux, flux_p, flux_m, u, u_p, u_m, E_p, E_m);
  
  for (int f = 0; f < 3; f++) {
    mrc_fld_destroy(flux_p[f]);
    mrc_fld_destroy(flux_m[f]);
    mrc_fld_destroy(u_delta[f]);
    mrc_fld_destroy(u_p[f]);
    mrc_fld_destroy(u_m[f]);
    mrc_fld_destroy(E_p[f]);
    mrc_fld_destroy(E_m[f]);
  }
  

}
