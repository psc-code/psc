
#ifndef PDE_MHD_CONVERT_C
#define PDE_MHD_CONVERT_C

#include "pde/pde_mhd_line.c"

#include "ggcm_mhd_private.h"

#define TINY_NUMBER 1.0e-20 // FIXME

// ----------------------------------------------------------------------
// convert_scons_from_prim

static inline void
convert_scons_from_prim(mrc_fld_data_t u[], mrc_fld_data_t w[])
{
  assert(s_gamma);
  
  u[RR ] = w[RR];
  u[RVX] = w[RR] * w[VX];
  u[RVY] = w[RR] * w[VY];
  u[RVZ] = w[RR] * w[VZ];
  u[UU ] = w[PP] * s_gamma_m1_inv +
    .5f * w[RR] * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ]));
  u[BX ] = w[BX];
  u[BY ] = w[BY];
  u[BZ ] = w[BZ];
}
      
// ----------------------------------------------------------------------
// mhd_pt_scons_from_prim

static inline void
mhd_pt_scons_from_prim(mrc_fld_data_t u[], mrc_fld_data_t w[])
{
  u[RR ] = w[RR];
  u[RVX] = w[RR] * w[VX];
  u[RVY] = w[RR] * w[VY];
  u[RVZ] = w[RR] * w[VZ];
  u[UU ] = w[PP] * s_gamma_m1_inv +
    .5 * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ])) * w[RR];
}

// ----------------------------------------------------------------------
// convert_prim_from_scons

static inline void
convert_prim_from_scons(mrc_fld_data_t w[], mrc_fld_data_t u[])
{
  assert(s_gamma);
  
  w[RR] = u[RR];
  mrc_fld_data_t rri = 1.f / u[RR];
  w[VX] = rri * u[RVX];
  w[VY] = rri * u[RVY];
  w[VZ] = rri * u[RVZ];
  mrc_fld_data_t rvv = (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri;
  w[PP] = s_gamma_m1 * (u[UU] - .5f * rvv);
  w[BX] = u[BX];
  w[BY] = u[BY];
  w[BZ] = u[BZ];
}
      
// ----------------------------------------------------------------------
// mhd_pt_prim_from_scons

static inline void
mhd_pt_prim_from_scons(mrc_fld_data_t w[], mrc_fld_data_t u[])
{
  w[RR] = u[RR];
  mrc_fld_data_t rri = 1.f / u[RR];
  w[VX] = rri * u[RVX];
  w[VY] = rri * u[RVY];
  w[VZ] = rri * u[RVZ];
  mrc_fld_data_t rvv = (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri;
  w[PP] = s_gamma_m1 * (u[UU] - .5f * rvv);
  w[PP] = mrc_fld_max(w[PP], TINY_NUMBER);
}

// ----------------------------------------------------------------------
// convert_fcons_from_prim

static inline void
convert_fcons_from_prim(mrc_fld_data_t u[], mrc_fld_data_t w[])
{
  assert(s_gamma);

  u[RR ] = w[RR];
  u[RVX] = w[RR] * w[VX];
  u[RVY] = w[RR] * w[VY];
  u[RVZ] = w[RR] * w[VZ];
  u[EE ] = w[PP] * s_gamma_m1_inv +
    .5f * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ])) * w[RR] +
    .5f * (sqr(w[BX]) + sqr(w[BY]) + sqr(w[BZ])) * s_mu0_inv;
  u[BX ] = w[BX];
  u[BY ] = w[BY];
  u[BZ ] = w[BZ];
  if (s_n_state == 9) {
    u[PSI] = 0.f;
  }
}

// ----------------------------------------------------------------------
// mhd_pt_fcons_from_prim

static inline void
mhd_pt_fcons_from_prim(mrc_fld_data_t u[], mrc_fld_data_t w[])
{
  u[RR ] = w[RR];
  u[RVX] = w[RR] * w[VX];
  u[RVY] = w[RR] * w[VY];
  u[RVZ] = w[RR] * w[VZ];
  u[EE ] = w[PP] * s_gamma_m1_inv +
    .5f * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ])) * w[RR] +
    .5f * (sqr(w[BX]) + sqr(w[BY]) + sqr(w[BZ])) * s_mu0_inv;
  u[BX ] = w[BX];
  u[BY ] = w[BY];
  u[BZ ] = w[BZ];
  if (s_opt_divb == OPT_DIVB_GLM) {
    u[PSI] = w[PSI];
  }
}

// ----------------------------------------------------------------------
// convert_prim_from_fcons

static inline void
convert_prim_from_fcons(mrc_fld_data_t w[], mrc_fld_data_t u[])
{
  assert(s_gamma);
  
  mrc_fld_data_t rri = 1.f / u[RR];
  w[RR] = u[RR];
  w[VX] = u[RVX] * rri;
  w[VY] = u[RVY] * rri;
  w[VZ] = u[RVZ] * rri;
  w[PP] = s_gamma_m1 * (u[EE] 
			- .5f * (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri
			- .5f * (sqr(u[BX] ) + sqr(u[BY] ) + sqr(u[BZ] )) * s_mu0_inv);
  w[BX] = u[BX];
  w[BY] = u[BY];
  w[BZ] = u[BZ];
}

// ----------------------------------------------------------------------
// mhd_pt_prim_from_fcons

static inline void
mhd_pt_prim_from_fcons(mrc_fld_data_t w[], mrc_fld_data_t u[])
{
  mrc_fld_data_t rri = 1.f / u[RR];
  w[RR] = u[RR];
  w[VX] = u[RVX] * rri;
  w[VY] = u[RVY] * rri;
  w[VZ] = u[RVZ] * rri;
  w[PP] = s_gamma_m1 * (u[EE] 
			- .5f * (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri
			- .5f * (sqr(u[BX] ) + sqr(u[BY] ) + sqr(u[BZ] )) * s_mu0_inv);
  w[PP] = mrc_fld_max(w[PP], TINY_NUMBER);
  w[BX] = u[BX];
  w[BY] = u[BY];
  w[BZ] = u[BZ];
  if (s_opt_divb == OPT_DIVB_GLM) {
    w[PSI] = u[PSI];
  }
}

// ----------------------------------------------------------------------
// convert_gkeyll_from_prim

static inline void
convert_gkeyll_from_prim(mrc_fld_data_t state[], mrc_fld_data_t prim[8])
{
  for (int s = 0; s < s_gk_nr_fluids; s++) {
    mrc_fld_data_t *state_sp = state + s_gk_idx[s];
    mrc_fld_data_t rrs = prim[RR] * s_gk_mass_ratios[s];
    state_sp[G5M_RRS ] = rrs;
    state_sp[G5M_RVXS] = rrs * prim[VX];
    state_sp[G5M_RVYS] = rrs * prim[VY];
    state_sp[G5M_RVZS] = rrs * prim[VZ];
    state_sp[G5M_UUS ] = (prim[PP] * s_gk_pressure_ratios[s]) * s_gamma_m1_inv
      + .5f * rrs * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
  }

  mrc_fld_data_t *state_em = state + s_gk_idx_em;
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
// convert_prim_from_gkeyll

static inline void
convert_prim_from_gkeyll(mrc_fld_data_t prim[8], mrc_fld_data_t state[])
{
  prim[RR] = 0.f;
  prim[VX] = prim[VY] = prim[VZ] = 0.f;
  prim[PP] = 0.f;
  for (int sp = 0; sp < s_gk_nr_fluids; sp++) {
    mrc_fld_data_t *state_sp = state + s_gk_idx[sp];
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
convert_state_from_prim(mrc_fld_data_t state[], mrc_fld_data_t prim[])
{
  if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS) {
    convert_scons_from_prim(state, prim);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
    convert_fcons_from_prim(state, prim);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL) {
    convert_gkeyll_from_prim(state, prim);
  } else {
    assert(0);
  }
}
      
// ----------------------------------------------------------------------
// convert_prim_from_state
      
static inline void
convert_prim_from_state(mrc_fld_data_t prim[], mrc_fld_data_t state[])
{
  if (MT_FORMULATION(MT) == MT_FORMULATION_SCONS) {
    convert_prim_from_scons(prim, state);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_FCONS) {
    convert_prim_from_fcons(prim, state);
  } else if (MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL) {
    convert_prim_from_gkeyll(prim, state);
  } else {
    assert(0);
  }
}

#endif

// ----------------------------------------------------------------------
// mhd_prim_from_fcons

static void
mhd_prim_from_fcons(fld1d_state_t W, fld1d_state_t U, int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    mhd_pt_prim_from_fcons(&F1S(W, 0, i), &F1S(U, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_fcons_from_prim

static void _mrc_unused
mhd_fcons_from_prim(fld1d_state_t U, fld1d_state_t W, int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    mhd_pt_fcons_from_prim(&F1S(U, 0, i), &F1S(W, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_prim_from_scons

static void
mhd_prim_from_scons(fld1d_state_t W, fld1d_state_t U, int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    mhd_pt_prim_from_scons(&F1S(W, 0, i), &F1S(U, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_scons_from_prim

static void _mrc_unused
mhd_scons_from_prim(fld1d_state_t U, fld1d_state_t W, int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    mhd_pt_scons_from_prim(&F1S(U, 0, i), &F1S(W, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_prim_from_cons

static void _mrc_unused
mhd_prim_from_cons(fld1d_state_t W, fld1d_state_t U, int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_prim_from_fcons(W, U, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_prim_from_scons(W, U, ib, ie);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_cons_from_prim

static void _mrc_unused
mhd_cons_from_prim(fld1d_state_t U, fld1d_state_t W, int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_fcons_from_prim(U, W, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_scons_from_prim(U, W, ib, ie);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_prim_from_cons

static void _mrc_unused
patch_prim_from_cons(fld3d_t p_W, fld3d_t p_U, int sw)
{
  static fld1d_state_t l_U, l_W;
  if (!fld1d_state_is_setup(l_U)) {
    fld1d_state_setup(&l_U);
    fld1d_state_setup(&l_W);
  }

  int dir = 0;
  pde_for_each_line(dir, j, k, sw) {
    int ib = -sw, ie = s_ldims[0] + sw;
    mhd_line_get_state(l_U, p_U, j, k, dir, ib, ie);
    mhd_prim_from_cons(l_W, l_U, ib, ie);
    mhd_line_put_state(l_W, p_W, j, k, dir, ib, ie);
  }
}

// ----------------------------------------------------------------------
// patch_prim_from_cons_v2
//
// FIXME: this is here for reference, and should eventually go away
// however, it may also have better performance
// and it gives exactly the same answer as Fortran

static void _mrc_unused
patch_prim_from_cons_v2(fld3d_t p_W, fld3d_t p_U, int sw)
{
  fld3d_foreach(i,j,k, sw, sw) {
    F3S(p_W, RR, i,j,k) = F3S(p_U, RR, i,j,k);
    mrc_fld_data_t rri  = 1.f / F3S(p_U, RR, i,j,k);
    F3S(p_W, VX, i,j,k) = rri * F3S(p_U, RVX, i,j,k);
    F3S(p_W, VY, i,j,k) = rri * F3S(p_U, RVY, i,j,k);
    F3S(p_W, VZ, i,j,k) = rri * F3S(p_U, RVZ, i,j,k);
    mrc_fld_data_t rvv =
      F3S(p_W, VX, i,j,k) * F3S(p_U, RVX, i,j,k) +
      F3S(p_W, VY, i,j,k) * F3S(p_U, RVY, i,j,k) +
      F3S(p_W, VZ, i,j,k) * F3S(p_U, RVZ, i,j,k);
    F3S(p_W, PP, i,j,k) = s_gamma_m1 * (F3S(p_U, UU, i,j,k) - .5f * rvv);
  } fld3d_foreach_end;
}

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

static inline void
convert_get_cc_state_from_3d(mrc_fld_data_t state[], struct mrc_fld *fld,
			     int i, int j, int k, int p)
{
#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  convert_get_state_from_3d(state, fld, i,j,k, p);
#else
  for (int m = 0; m < 5; m++) {
    state[m] = M3(fld, m, i,j,k, p);
  }
#if MT_BGRID(MT) == MT_BGRID_FC
  state[BX] = .5f * (BX_(fld, i,j,k, p) + BX_(fld, i+di,j   ,k   , p));
  state[BY] = .5f * (BY_(fld, i,j,k, p) + BY_(fld, i   ,j+dj,k   , p));
  state[BZ] = .5f * (BZ_(fld, i,j,k, p) + BZ_(fld, i   ,j   ,k+dk, p));
#elif MT_BGRID(MT) == MT_BGRID_FC_GGCM
  state[BX] = .5f * (BX_(fld, i,j,k, p) + BX_(fld, i-di,j   ,k   , p));
  state[BY] = .5f * (BY_(fld, i,j,k, p) + BY_(fld, i   ,j-dj,k   , p));
  state[BZ] = .5f * (BZ_(fld, i,j,k, p) + BZ_(fld, i   ,j   ,k-dk, p));
#elif MT_BGRID(MT) == MT_BGRID_CC
  state[BX] = BX_(fld, i,j,k, p);
  state[BY] = BY_(fld, i,j,k, p);
  state[BZ] = BZ_(fld, i,j,k, p);
#endif
#endif
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
  for (int sp = 0; sp < s_gk_nr_fluids; sp++) {
    for (int m = s_gk_idx[sp]; m < s_gk_idx[sp] + G5M_NRS; m++) {
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
