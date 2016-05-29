
#include "ggcm_mhd_private.h"

#define TINY_NUMBER 1.0e-20 // FIXME

// ----------------------------------------------------------------------
// mhd_prim_from_fc

static void _mrc_unused
mhd_prim_from_fc(struct ggcm_mhd *mhd, struct mrc_fld *W_cc, struct mrc_fld *U_cc,
		 int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = mhd->par.gamm - 1.;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *W = &F1(W_cc, 0, i), *U = &F1(U_cc, 0, i);

    mrc_fld_data_t rri = 1. / U[RR];
    W[RR] = U[RR];
    W[VX] = U[RVX] * rri;
    W[VY] = U[RVY] * rri;
    W[VZ] = U[RVZ] * rri;
    W[PP] = gamma_minus_1 * (U[EE] 
			     - .5 * (sqr(U[RVX]) + sqr(U[RVY]) + sqr(U[RVZ])) * rri
			     - .5 * (sqr(U[BX]) + sqr(U[BY]) + sqr(U[BZ])));
    W[PP] = fmax(W[PP], TINY_NUMBER);
    W[BX] = U[BX];
    W[BY] = U[BY];
    W[BZ] = U[BZ];
  }
}

// ----------------------------------------------------------------------
// mhd_fc_from_prim

static void _mrc_unused
mhd_fc_from_prim(struct ggcm_mhd *mhd, struct mrc_fld *U_cc, struct mrc_fld *W_cc,
		    int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = mhd->par.gamm - 1.;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *U = &F1(U_cc, 0, i), *W = &F1(W_cc, 0, i);

    mrc_fld_data_t rr = W[RR];
    U[RR ] = rr;
    U[RVX] = rr * W[VX];
    U[RVY] = rr * W[VY];
    U[RVZ] = rr * W[VZ];
    U[EE ] = 
      W[PP] / gamma_minus_1 +
      + .5 * (sqr(W[VX]) + sqr(W[VY]) + sqr(W[VZ])) * rr
      + .5 * (sqr(W[BX]) + sqr(W[BY]) + sqr(W[BZ]));
    U[BX ] = W[BX];
    U[BY ] = W[BY];
    U[BZ ] = W[BZ];
  }
}

