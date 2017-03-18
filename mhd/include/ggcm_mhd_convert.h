
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
// ----------------------------------------------------------------------	
// convert_primitive_5m_point_comove

static inline void
convert_primitive_5m_point_comove(mrc_fld_data_t state[], mrc_fld_data_t vals[], int nr_fluids, int nr_moments,
    float mass[], float charge[], float pressure_ratios[], float gamm)
{
  mrc_fld_data_t mass_ratios[nr_fluids];
  mrc_fld_data_t mass_total = 0.;
  int idx[nr_fluids];

  for (int s = 0; s < nr_fluids; s++) {
    mass_total += mass[s];
    idx[s] = s * nr_moments;
  }
  for (int s = 0; s < nr_fluids; s++)
    mass_ratios[s] = mass[s] / mass_total;

  int idx_em = nr_fluids * nr_moments;

  mrc_fld_data_t rr = vals[RR];
  mrc_fld_data_t vx = vals[VX];
  mrc_fld_data_t vy = vals[VY];
  mrc_fld_data_t vz = vals[VZ];
  mrc_fld_data_t pp = vals[PP];
  mrc_fld_data_t bx = vals[BX];
  mrc_fld_data_t by = vals[BY];
  mrc_fld_data_t bz = vals[BZ];

  for (int s = 0; s < nr_fluids; s++) {
    state[idx[s] + G5M_RRS ] = rr * mass_ratios[s];
    state[idx[s] + G5M_RVXS] = rr * mass_ratios[s] * vx;
    state[idx[s] + G5M_RVYS] = rr * mass_ratios[s] * vy;
    state[idx[s] + G5M_RVZS] = rr * mass_ratios[s] * vz;
    state[idx[s] + G5M_UUS ] = pp * pressure_ratios[s] / (gamm - 1.)
      + .5 * (sqr(state[idx[s] + G5M_RVXS])
            + sqr(state[idx[s] + G5M_RVYS])
            + sqr(state[idx[s] + G5M_RVZS])) / state[idx[s] + G5M_RRS];
  }

  state[idx_em + GK_EX] = - vy * bz + vz * by;
  state[idx_em + GK_EY] = - vz * bx + vx * bz;
  state[idx_em + GK_EZ] = - vx * by + vy * bx;

  state[idx_em + GK_BX] = bx;
  state[idx_em + GK_BY] = by;
  state[idx_em + GK_BZ] = bz;

  state[idx_em + GK_PHI] = 0.;
  state[idx_em + GK_PSI] = 0.;
}

// ----------------------------------------------------------------------
// convert_state_from_prim_gkeyll

static inline void
convert_state_from_prim_gkeyll(mrc_fld_data_t state[], mrc_fld_data_t prim[8])
{
  // FIXME: partitioning of mhd quantities between species is probably too rough
  
  for (int sp = 0; sp < cvt_gk_nr_fluids; sp++) {
    mrc_fld_data_t *state_sp = state + cvt_gk_idx[sp];
    mrc_fld_data_t rrs = prim[RR] * cvt_gk_mass_ratios[sp];
    state_sp[G5M_RRS ] = rrs;
    state_sp[G5M_RVXS] = rrs * prim[VX];
    state_sp[G5M_RVYS] = rrs * prim[VY];
    state_sp[G5M_RVZS] = rrs * prim[VZ];
    state_sp[G5M_UUS ] = (prim[PP] * cvt_gk_pressure_ratios[sp]) * cvt_gamma_m1_inv
      + .5f * rrs * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  }
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

#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
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
#endif

// ----------------------------------------------------------------------
// convert_put_fluid_state_to_3d

static inline void
convert_put_fluid_state_to_3d(mrc_fld_data_t state[], struct mrc_fld *f,
			      int i, int j, int k, int p)
{
#if MT_FORMULATION(MT) == MT_FORMULATION_SCONS || MT_FORMULATION(MT) == MT_FORMULATION_FCONS
  convert_put_fluid_state_to_3d_mhd(state, f, i,j,k, p);
#elif MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  convert_put_fluid_state_to_3d_gkeyll(state, f, i,j,k, p);
#else
  assert(0);
#endif
}

#endif
