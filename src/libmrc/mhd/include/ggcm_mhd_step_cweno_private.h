
#ifndef GGCM_MHD_STEP_CWENO_PRIVATE_H
#define GGCM_MHD_STEP_CWENO_PRIVATE_H

#include "ggcm_mhd_step.h"
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_defs.h"

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_io.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include <mrc_fld.h>

#define FLUX(flux, f, m, ix,iy,iz) MRC_F3((flux)[f], m, ix,iy,iz)
#define sign(x) (( x > 0 ) - ( x < 0 ))

// Toggle debugging output/tests
//#define DEBUG

// Define limiter: [0] none [1] van Leer, [2] minmod [3] moncen [4] genminmod [5] Van Albada
#define LMTR 2

// KNP[0] or KT[1]?
#define KT 0

// include whistler speed?
#define incws 0

// toggle between semi-conservative and conservative form
#define SEMICONSV 1

// reduction factor for JxB and pressure terms
#define RFACT 1.0
//#define RFACT MRC_F3(mhd->fld, _ZMASK, ix,iy,iz)

//density floor
#define RMIN 6.000000e-04

//SW BND toggle
#define SWBND 0

//CWENO REC for fluid variables toggle
#define CWENOREC 1

enum {
  // reuse B in the _fluxes_ (only) to store E field
  _EX = BX,
  _EY,
  _EZ,

  _JX = 11,
  _JY,
  _JZ,

  __NR_FLDS,
};

// function prototypes:


// ----------------------------------------------------------------------
// ggcm_mhd_get_fields

struct mrc_fld * ggcm_mhd_get_fields(struct ggcm_mhd *mhd, const char *name, int nr_comp);

// ----------------------------------------------------------------------
// fill_ghost_fld
//
// This fills ghost cells for fld objects that have been duplicated in the
// time-stepping routine (c.f. mrc_ts_rk2.c). Without this, zero-values at
// boundaries will give inf and nan values for non-periodic boudnaries. (obsolete)

void ggcm_mhd_fill_ghost_fld(struct ggcm_mhd *mhd, struct mrc_fld *_fld);

// ----------------------------------------------------------------------
// calc_u_delta
//
// calculate reconstruction to cell interfaces

void calc_u_delta(struct mrc_fld *_u_delta[3], struct mrc_fld *_u, struct ggcm_mhd *mhd);

// ----------------------------------------------------------------------
// calc_semiconsv_rhs
//
// calculates rhs for semi-conservative mhd

void calc_semiconsv_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3]);

// ----------------------------------------------------------------------
// calc_neg_divg
//
// calculates negative divergence

void calc_neg_divg(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3]);

// ----------------------------------------------------------------------
// calc_fluxes_per_face
//
// this calculates fluxes on the face i using reconstructed variables (fld) that are
// given on the respective face
// (Ziegler 2004 section 3.1)

void calc_fluxes_per_face(struct mrc_fld **_flux, struct ggcm_mhd *mhd, struct mrc_fld *_fld, int i);


// ----------------------------------------------------------------------
// calc_u_pm
// Calculate point values of variables at cell surface centers using
// linear reconstrcution with TVD slope limiters
// (Ziegler 2004 section 3.2)
void calc_u_pm(struct ggcm_mhd *mhd, struct mrc_fld *_u_p[3], struct mrc_fld *_u_m[3],
	       struct mrc_fld *_E_p[3], struct mrc_fld *_E_m[3],
	       struct mrc_fld *_u, struct mrc_fld *_u_delta[3]);


// ----------------------------------------------------------------------
// calc_KNP_fluxes
// A. Kurganov, S. Noelle, G. Petrova, SIAM J. Sci. Comput. 23 (2001) 707.
// (Ziegler 2004 section 3.1)

void calc_KNP_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
		     struct mrc_fld *_flux_p[3], struct mrc_fld *_flux_m[3],
		     struct mrc_fld *_u,
		     struct mrc_fld *_u_p[3], struct mrc_fld *_u_m[3],
		     struct mrc_fld *_E_p[3], struct mrc_fld *_E_m[3]);

// ----------------------------------------------------------------------
// calc_ct_rhs
//
// calculate induction equation rhs with constrained transport method
// (Ziegler 2004 section 3.3)
void calc_ct_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3]);


// ----------------------------------------------------------------------
// calc_cweno_fluxes
//
// calculates CWENO fluxes on faces in flux_E, from the original state
// vector u (which is cell centered / on the Yee grid)
// flux[0-4] are the fluid vars
// flux[5-7] are E-field "fluxes" (not B-field!)

void calc_cweno_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
		       struct mrc_fld *_u);

void calc_u_cweno(struct ggcm_mhd *mhd, struct mrc_fld *u_p[3], struct mrc_fld *u_m[3],
                  struct mrc_fld *E_p[3], struct mrc_fld *E_m[3],
                  struct mrc_fld *u, struct mrc_fld *u_delta[3]);

#endif
