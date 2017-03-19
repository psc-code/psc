
#ifndef GGCM_MHD_CONVERT_H
#define GGCM_MHD_CONVERT_H

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_gkeyll.h"

#include "pde/pde_defs.h"

#define OPT_FLD1D OPT_FLD1D_C_ARRAY

#include "pde/pde_mhd_setup.c"

static inline void
convert_state_from_prim_scons(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  assert(s_gamma);
  
  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[UU ] = prim[PP] * s_gamma_m1_inv +
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  state[BX ] = prim[BX];
  state[BY ] = prim[BY];
  state[BZ ] = prim[BZ];
}
      
static inline void
convert_prim_from_state_scons(mrc_fld_data_t prim[8], mrc_fld_data_t state[8])
{
  assert(s_gamma);
  
  prim[RR] = state[RR];
  mrc_fld_data_t rri = 1.f / state[RR];
  prim[VX] = rri * state[RVX];
  prim[VY] = rri * state[RVY];
  prim[VZ] = rri * state[RVZ];
  prim[PP] = s_gamma_m1 * (state[UU] 
			   - .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ])));
  prim[BX] = state[BX];
  prim[BY] = state[BY];
  prim[BZ] = state[BZ];
}
      
static inline void
convert_state_from_prim_fcons(mrc_fld_data_t state[], mrc_fld_data_t prim[8])
{
  assert(s_gamma);

  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[EE ] = prim[PP] * s_gamma_m1_inv
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]))
    + .5f * (sqr(prim[BX]) + sqr(prim[BY]) + sqr(prim[BZ]));
  state[BX ] = prim[BX];
  state[BY ] = prim[BY];
  state[BZ ] = prim[BZ];
  if (s_n_state == 9) {
    state[PSI] = 0.f;
  }
}

static inline void
convert_prim_from_state_fcons(mrc_fld_data_t prim[8], mrc_fld_data_t state[])
{
  assert(s_gamma);
  
  prim[RR] = state[RR];
  mrc_fld_data_t rri = 1.f / state[RR];
  prim[VX] = rri * state[RVX];
  prim[VY] = rri * state[RVY];
  prim[VZ] = rri * state[RVZ];
  prim[PP] = s_gamma_m1 * (state[UU] 
			   - .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]))
			   - .5f * (sqr(state[BX]) + sqr(state[BY]) + sqr(state[BZ])));
  prim[BX] = state[BX];
  prim[BY] = state[BY];
  prim[BZ] = state[BZ];
}

// ----------------------------------------------------------------------
// convert_state_from_prim_gkeyll

static inline void
convert_state_from_prim_gkeyll(mrc_fld_data_t state[], mrc_fld_data_t prim[8])
{
  for (int s = 0; s < cvt_gk_nr_fluids; s++) {
    mrc_fld_data_t *state_sp = state + cvt_gk_idx[s];
    mrc_fld_data_t rrs = prim[RR] * cvt_gk_mass_ratios[s];
    state_sp[G5M_RRS ] = rrs;
    state_sp[G5M_RVXS] = rrs * prim[VX];
    state_sp[G5M_RVYS] = rrs * prim[VY];
    state_sp[G5M_RVZS] = rrs * prim[VZ];
    state_sp[G5M_UUS ] = (prim[PP] * cvt_gk_pressure_ratios[s]) * s_gamma_m1_inv
      + .5f * rrs * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  }

  mrc_fld_data_t *state_em = state + cvt_gk_idx_em;
  state_em[GK_EX] = - prim[VY] * prim[BZ] + prim[VZ] * prim[BY];
  state_em[GK_EY] = - prim[VZ] * prim[BX] + prim[VX] * prim[BZ];
  state_em[GK_EZ] = - prim[VX] * prim[BY] + prim[VY] * prim[BX];

  state_em[GK_BX] = prim[BX];
  state_em[GK_BY] = prim[BY];
  state_em[GK_BZ] = prim[BZ];

  state_em[GK_PHI] = 0.;
  state_em[GK_PSI] = 0.;
}

// ----------------------------------------------------------------------
// convert_prim_from_state_gkeyll

static inline void
convert_prim_from_state_gkeyll(mrc_fld_data_t prim[8], mrc_fld_data_t state[])
{
  prim[RR] = 0.f;
  prim[VX] = prim[VY] = prim[VZ] = 0.f;
  prim[PP] = 0.f;
  for (int sp = 0; sp < cvt_gk_nr_fluids; sp++) {
    mrc_fld_data_t *state_sp = state + cvt_gk_idx[sp];
    mrc_fld_data_t rrs = state_sp[G5M_RRS];
    prim[RR] += rrs;
    prim[VX] += state_sp[G5M_RVXS];
    prim[VY] += state_sp[G5M_RVYS];
    prim[VZ] += state_sp[G5M_RVZS];
    mrc_fld_data_t rvvs = (sqr(state_sp[G5M_RVXS]) +
			   sqr(state_sp[G5M_RVYS]) +
			   sqr(state_sp[G5M_RVZS])) / rrs;
    prim[PP] += s_gamma_m1 * (state_sp[G5M_UUS] - .5f * rvvs);
  }
  prim[VX] /= prim[RR];
  prim[VY] /= prim[RR];
  prim[VZ] /= prim[RR];
}

#ifdef MT

// ----------------------------------------------------------------------
// convert_state_from_prim
      
static inline void
convert_state_from_prim(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS) {
    convert_state_from_prim_scons(state, prim);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
    convert_state_from_prim_fcons(state, prim);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL) {
    convert_state_from_prim_gkeyll(state, prim);
  } else {
    assert(0);
  }
}
      
// ----------------------------------------------------------------------
// convert_prim_from_state
      
static inline void
convert_prim_from_state(mrc_fld_data_t prim[8], mrc_fld_data_t state[8])
{
  if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS) {
    convert_prim_from_state_scons(prim, state);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
    convert_prim_from_state_fcons(prim, state);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL) {
    convert_prim_from_state_gkeyll(prim, state);
  } else {
    assert(0);
  }
}

#endif

// ----------------------------------------------------------------------
// convert_get_state_from_3d
      
static inline void
convert_get_state_from_3d(mrc_fld_data_t state[], struct mrc_fld *f,
			  int i, int j, int k, int p)
{
  for (int m = 0; m < s_n_state; m++) {
    state[m] = M3(f, m, i,j,k, p);
  }
}

// ----------------------------------------------------------------------
// convert_put_fluid_state_to_3d_mhd

static inline void
convert_put_fluid_state_to_3d_mhd(mrc_fld_data_t state[], struct mrc_fld *f,
				  int i, int j, int k, int p)
{
  for (int m = 0; m < 5; m ++) {
    M3(f, m, i,j,k, p) = state[m];
  }
}

// ----------------------------------------------------------------------
// convert_put_fluid_state_to_3d_gkeyll

static inline void
convert_put_fluid_state_to_3d_gkeyll(mrc_fld_data_t state[], struct mrc_fld *f,
				     int i, int j, int k, int p)
{
  for (int sp = 0; sp < cvt_gk_nr_fluids; sp++) {
    for (int m = cvt_gk_idx[sp]; m < cvt_gk_idx[sp] + G5M_NRS; m++) {
      M3(f, m, i,j,k, p) = state[m];
    }
  }
}

#ifdef MT

// ----------------------------------------------------------------------
// convert_put_fluid_state_to_3d

static inline void
convert_put_fluid_state_to_3d(mrc_fld_data_t state[], struct mrc_fld *f,
			      int i, int j, int k, int p)
{
  if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS ||
      MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
    convert_put_fluid_state_to_3d_mhd(state, f, i,j,k, p);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL) {
    convert_put_fluid_state_to_3d_gkeyll(state, f, i,j,k, p);
  } else {
    assert(0);
  }
}

#endif

// ----------------------------------------------------------------------
// convert_put_state_to_3d

static inline void
convert_put_state_to_3d(mrc_fld_data_t state[], struct mrc_fld *f,
			int i, int j, int k, int p)
{
  for (int m = 0; m < s_n_state; m++)
    M3(f, m, i,j,k, p) = state[m];
}

#endif
