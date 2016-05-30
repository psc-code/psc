
#include "ggcm_mhd_private.h"

#define TINY_NUMBER 1.0e-20 // FIXME

// ----------------------------------------------------------------------
// mhd_prim_from_fcons

static void
mhd_prim_from_fcons(fld1d_state_t W, fld1d_state_t U, int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = s_gamma - 1.f;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *w = &F1S(W, 0, i), *u = &F1S(U, 0, i);

    mrc_fld_data_t rri = 1.f / u[RR];
    w[RR] = u[RR];
    w[VX] = u[RVX] * rri;
    w[VY] = u[RVY] * rri;
    w[VZ] = u[RVZ] * rri;
    w[PP] = gamma_minus_1 * (u[EE] 
			     - .5 * (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri
			     - .5 * (sqr(u[BX]) + sqr(u[BY]) + sqr(u[BZ])));
    w[PP] = fmax(w[PP], TINY_NUMBER);
    w[BX] = u[BX];
    w[BY] = u[BY];
    w[BZ] = u[BZ];
  }
}

// ----------------------------------------------------------------------
// mhd_fc_from_prim

static void _mrc_unused
mhd_fc_from_prim(fld1d_state_t U, fld1d_state_t W, int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = s_gamma - 1.f;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *u = &F1S(U, 0, i), *w = &F1S(W, 0, i);

    mrc_fld_data_t rr = w[RR];
    u[RR ] = rr;
    u[RVX] = rr * w[VX];
    u[RVY] = rr * w[VY];
    u[RVZ] = rr * w[VZ];
    u[EE ] = 
      w[PP] / gamma_minus_1 +
      + .5 * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ])) * rr
      + .5 * (sqr(w[BX]) + sqr(w[BY]) + sqr(w[BZ]));
    u[BX ] = w[BX];
    u[BY ] = w[BY];
    u[BZ ] = w[BZ];
  }
}

// ----------------------------------------------------------------------
// mhd_prim_from_scons

static void
mhd_prim_from_scons(fld1d_state_t W, fld1d_state_t U, int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = s_gamma - 1.f;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *w = &F1S(W, 0, i), *u = &F1S(U, 0, i);

    mrc_fld_data_t rri = 1.f / u[RR];
    w[RR] = u[RR];
    w[VX] = rri * u[RVX];
    w[VY] = rri * u[RVY];
    w[VZ] = rri * u[RVZ];
    mrc_fld_data_t rvv = (sqr(u[RVX]) + sqr(u[RVY]) + sqr(u[RVZ])) * rri;
    w[PP] = gamma_minus_1 * (u[UU] - .5 * rvv);
    w[PP] = fmax(w[PP], TINY_NUMBER);
  }
}

// ----------------------------------------------------------------------
// mhd_sc_from_prim

static void _mrc_unused
mhd_sc_from_prim(fld1d_state_t U, fld1d_state_t W, int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = s_gamma - 1.f;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *u = &F1S(U, 0, i), *w = &F1S(W, 0, i);

    mrc_fld_data_t rr = w[RR];
    u[RR ] = rr;
    u[RVX] = rr * w[VX];
    u[RVY] = rr * w[VY];
    u[RVZ] = rr * w[VZ];
    u[UU ] = 
      w[PP] / gamma_minus_1 +
      + .5 * (sqr(w[VX]) + sqr(w[VY]) + sqr(w[VZ])) * rr;
  }
}

// ----------------------------------------------------------------------
// mhd_prim_from_cons

static void _mrc_unused
mhd_prim_from_cons(fld1d_state_t W, fld1d_state_t U, int ldim, int l, int r)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_prim_from_fcons(W, U, ldim, l, r);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_prim_from_scons(W, U, ldim, l, r);
  } else {
    assert(0);
  }
}
