
#ifndef GGCM_MHD_CONVERT_H
#define GGCM_MHD_CONVERT_H

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"

static mrc_fld_data_t cvt_gamma;
static mrc_fld_data_t cvt_gamma_m1;
static mrc_fld_data_t cvt_gamma_m1_inv;
static int cvt_n_state;
#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
static int cvt_gk_nr_fluids;
static int *cvt_gk_idx;
static int cvt_gk_idx_em;
static float *cvt_gk_mass_ratios;
static float *cvt_gk_pressure_ratios;
#endif

static inline void ggcm_mhd_convert_setup(struct ggcm_mhd *mhd)
{
  cvt_gamma = mhd->par.gamm;
  cvt_gamma_m1 = cvt_gamma - 1.f;
  cvt_gamma_m1_inv = 1.f / cvt_gamma_m1;
  cvt_n_state = mrc_fld_nr_comps(mhd->fld);
  // FIXME, hacky as usual, to deal with the legacy all-in-one big array
  if (cvt_n_state == _NR_FLDS) {
    cvt_n_state = 8;
  }

#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  cvt_gk_nr_fluids = mhd->par.gk_nr_fluids;
  cvt_gk_idx = mhd->par.gk_idx;
  cvt_gk_idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);
  cvt_gk_mass_ratios = mhd->par.gk_mass_ratios;
  cvt_gk_pressure_ratios = mhd->par.gk_pressure_ratios.vals;
#endif
}

static inline void
convert_state_from_prim_scons(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  assert(cvt_gamma);
  
  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[UU ] = prim[PP] * cvt_gamma_m1_inv +
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  state[BX ] = prim[BX];
  state[BY ] = prim[BY];
  state[BZ ] = prim[BZ];
}
      
static inline void
convert_prim_from_state_scons(mrc_fld_data_t prim[8], mrc_fld_data_t state[8])
{
  assert(cvt_gamma);
  
  prim[RR] = state[RR];
  mrc_fld_data_t rri = 1.f / state[RR];
  prim[VX] = rri * state[RVX];
  prim[VY] = rri * state[RVY];
  prim[VZ] = rri * state[RVZ];
  prim[PP] = cvt_gamma_m1 * (state[UU] 
			     - .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ])));
  prim[BX] = state[BX];
  prim[BY] = state[BY];
  prim[BZ] = state[BZ];
}
      
static inline void
convert_state_from_prim_fcons(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  assert(cvt_gamma);

  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[EE ] = prim[PP] * cvt_gamma_m1_inv
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]))
    + .5f * (sqr(prim[BX]) + sqr(prim[BY]) + sqr(prim[BZ]));
  state[BX ] = prim[BX];
  state[BY ] = prim[BY];
  state[BZ ] = prim[BZ];
}

static inline void
convert_prim_from_state_fcons(mrc_fld_data_t prim[8], mrc_fld_data_t state[8])
{
  assert(cvt_gamma);
  
  prim[RR] = state[RR];
  mrc_fld_data_t rri = 1.f / state[RR];
  prim[VX] = rri * state[RVX];
  prim[VY] = rri * state[RVY];
  prim[VZ] = rri * state[RVZ];
  prim[PP] = cvt_gamma_m1 * (state[UU] 
			     - .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]))
			     - .5f * (sqr(state[BX]) + sqr(state[BY]) + sqr(state[BZ])));
  prim[BX] = state[BX];
  prim[BY] = state[BY];
  prim[BZ] = state[BZ];
}

#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
static inline void
convert_state_from_prim_gkeyll(mrc_fld_data_t state[], mrc_fld_data_t prim[8])
{
  // FIXME: partitioning of mhd quantities between species is probably too rough
  
  for (int sp = 0; sp < cvt_gk_nr_fluids; sp++) {
    mrc_fld_data_t *state_sp = state + cvt_gk_idx[sp];
    float rrs = prim[RR] * cvt_gk_mass_ratios[sp];
    state_sp[G5M_RRS ] = rrs;
    state_sp[G5M_RVXS] = rrs * prim[VX];
    state_sp[G5M_RVYS] = rrs * prim[VY];
    state_sp[G5M_RVZS] = rrs * prim[VZ];
    state_sp[G5M_UUS ] = (prim[PP] * cvt_gk_pressure_ratios[sp]) * cvt_gamma_m1_inv
      + .5f * rrs * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  }
}

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
    prim[PP] += cvt_gamma_m1 * (state_sp[G5M_UUS] - .5f * rvvs);
  }
  prim[VX] /= prim[RR];
  prim[VY] /= prim[RR];
  prim[VZ] /= prim[RR];
}
#endif

// ----------------------------------------------------------------------
// convert_state_from_prim
      
static inline void
convert_state_from_prim(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
#if MT_FORMULATION(MT) == MT_FORMULATION_SCONS
  convert_state_from_prim_scons(state, prim);
#elif MT_FORMULATION(MT) == MT_FORMULATION_FCONS
  convert_state_from_prim_fcons(state, prim);
#elif MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  convert_state_from_prim_gkeyll(state, prim);
#else
  assert(0);
#endif
}
      
// ----------------------------------------------------------------------
// convert_prim_from_state
      
static inline void
convert_prim_from_state(mrc_fld_data_t prim[8], mrc_fld_data_t state[8])
{
#if MT_FORMULATION(MT) == MT_FORMULATION_SCONS
  convert_prim_from_state_scons(prim, state);
#elif MT_FORMULATION(MT) == MT_FORMULATION_FCONS
  convert_prim_from_state_fcons(prim, state);
#elif MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  convert_prim_from_state_gkeyll(prim, state);
#else
  assert(0);
#endif
}

// ----------------------------------------------------------------------
// convert_get_state_from_3d
      
static inline void
convert_get_state_from_3d(mrc_fld_data_t state[], struct mrc_fld *f,
			  int i, int j, int k, int p)
{
  for (int m = 0; m < cvt_n_state; m++) {
    state[m] = M3(f, m, i,j,k, p);
  }
}

#endif
