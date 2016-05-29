
#include "ggcm_mhd_private.h"

#include <math.h>

#define TINY_NUMBER 1.0e-20 // FIXME

// ======================================================================
// fully-conservative

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

// ======================================================================
// semi-conservative

// ----------------------------------------------------------------------
// mhd_sc_from_prim_1d

static inline void
mhd_sc_from_prim_1d(struct ggcm_mhd *mhd, struct mrc_fld *U1d, struct mrc_fld *W1d,
		    int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = mhd->par.gamm - 1.;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *U = &F1(U1d, 0, i), *W = &F1(W1d, 0, i);

    U[RR ] = W[RR];
    U[RVX] = W[RR] * W[VX];
    U[RVY] = W[RR] * W[VY];
    U[RVZ] = W[RR] * W[VZ];
    U[UU ] = 
      W[PP] / gamma_minus_1 +
      + .5 * (sqr(W[VX]) + sqr(W[VY]) + sqr(W[VZ])) * W[RR];
  }
}

// ----------------------------------------------------------------------
// mhd_prim_from_sc

static void _mrc_unused
mhd_prim_from_sc(struct ggcm_mhd *mhd, struct mrc_fld *W_cc, struct mrc_fld *U_cc,
		 int ldim, int l, int r)
{
  mrc_fld_data_t gamma_minus_1 = mhd->par.gamm - 1.;

  for (int i = -l; i < ldim + r; i++) {
    mrc_fld_data_t *W = &F1(W_cc, 0, i), *U = &F1(U_cc, 0, i);

    mrc_fld_data_t rri = 1. / U[RR];
    W[RR] = U[RR];
    W[VX] = rri * U[RVX];
    W[VY] = rri * U[RVY];
    W[VZ] = rri * U[RVZ];
    mrc_fld_data_t rvv = (sqr(U[RVX]) + sqr(U[RVY]) + sqr(U[RVZ])) * rri;
    W[PP] = gamma_minus_1 * (U[UU] - .5 * rvv);
    W[PP] = fmax(W[PP], TINY_NUMBER);
  }
}

// ----------------------------------------------------------------------
// mhd_sc_from_prim

static void _mrc_unused
mhd_sc_from_prim(struct ggcm_mhd *mhd, struct mrc_fld *U_cc, struct mrc_fld *W_cc,
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
      + .5 * (sqr(W[VX]) + sqr(W[VY]) + sqr(W[VZ])) * rr;
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// pick_line_fc

static void _mrc_unused
pick_line_fc(struct mrc_fld *U_1d, struct mrc_fld *Bxi,
	     struct mrc_fld *U, struct mrc_fld *Bcc,
	     int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PICK_LINE(X, Y, Z, I, J, K)			       \
  do {							       \
    for (int i = - l; i < ldim + r; i++) {		       \
      F1(U_1d, RR , i) = M3(U, RR     , I,J,K, p);	       \
      F1(U_1d, RVX, i) = M3(U, RVX + X, I,J,K, p);	       \
      F1(U_1d, RVY, i) = M3(U, RVX + Y, I,J,K, p);	       \
      F1(U_1d, RVZ, i) = M3(U, RVX + Z, I,J,K, p);	       \
      F1(U_1d, EE , i) = M3(U, EE     , I,J,K, p);	       \
      F1(U_1d, BX , i) = M3(Bcc, X    , I,J,K, p);	       \
      F1(U_1d, BY , i) = M3(Bcc, Y    , I,J,K, p);	       \
      F1(U_1d, BZ , i) = M3(Bcc, Z    , I,J,K, p);	       \
      F1(Bxi, 0, i)    = M3(U, BX+X   , I,J,K, p);	       \
    }							       \
  } while (0)
  if (dir == 0) {
    PICK_LINE(0, 1, 2, i, j, k);
  } else if (dir == 1) {
    PICK_LINE(1, 2, 0, k, i, j);
  } else if (dir == 2) {
    PICK_LINE(2, 0, 1, j, k, i);
  } else {
    assert(0);
  }
#undef PICK_LINE
}

// ----------------------------------------------------------------------
// pick_line_sc

static void _mrc_unused
pick_line_sc(struct mrc_fld *x1, struct mrc_fld *x,
	     int ldim, int l, int r, int j, int k, int dim, int p)
{
#define PICK_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = -l; i < ldim + r; i++) {		\
      F1(x1, RR , i) = M3(x, RR   , I,J,K, p);		\
      F1(x1, RVX, i) = M3(x, RVX+X, I,J,K, p);		\
      F1(x1, RVY, i) = M3(x, RVX+Y, I,J,K, p);		\
      F1(x1, RVZ, i) = M3(x, RVX+Z, I,J,K, p);		\
      F1(x1, UU , i) = M3(x, UU   , I,J,K, p);		\
    }							\
  } while (0)

  if (dim == 0) {
    PICK_LINE(0,1,2, i,j,k);
  } else if (dim == 1) {
    PICK_LINE(1,2,0, k,i,j);
  } else if (dim == 2) {
    PICK_LINE(2,0,1, j,k,i);
  }
#undef PICK_LINE
}

// ----------------------------------------------------------------------
// put_line_fc

static void _mrc_unused
put_line_fc(struct mrc_fld *flux, struct mrc_fld *F_1d,
	    int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = -l; i < ldim + r; i++) {			\
      M3(flux, RR     , I,J,K, p) = F1(F_1d, RR , i);		\
      M3(flux, RVX + X, I,J,K, p) = F1(F_1d, RVX, i);		\
      M3(flux, RVX + Y, I,J,K, p) = F1(F_1d, RVY, i);		\
      M3(flux, RVX + Z, I,J,K, p) = F1(F_1d, RVZ, i);		\
      M3(flux, EE     , I,J,K, p) = F1(F_1d, EE , i);		\
      M3(flux, BX + Y , I,J,K, p) = F1(F_1d, BY , i);		\
      M3(flux, BX + Z , I,J,K, p) = F1(F_1d, BZ , i);		\
    }								\
} while (0)

  if (dir == 0) {
    PUT_LINE(0, 1, 2, i, j, k);
  } else if (dir == 1) {
    PUT_LINE(1, 2, 0, k, i, j);
  } else if (dir == 2) {
    PUT_LINE(2, 0, 1, j, k, i);
  } else {
    assert(0);
  }
#undef PUT_LINE
}

// ----------------------------------------------------------------------
// put_line_fc

static void _mrc_unused
put_line_fc_cc(struct mrc_fld *flux, struct mrc_fld *F_1d,
	       int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = -l; i < ldim + r; i++) {			\
      M3(flux, RR     , I,J,K, p) = F1(F_1d, RR , i);		\
      M3(flux, RVX + X, I,J,K, p) = F1(F_1d, RVX, i);		\
      M3(flux, RVX + Y, I,J,K, p) = F1(F_1d, RVY, i);		\
      M3(flux, RVX + Z, I,J,K, p) = F1(F_1d, RVZ, i);		\
      M3(flux, EE     , I,J,K, p) = F1(F_1d, EE , i);		\
      M3(flux, BX + X , I,J,K, p) = F1(F_1d, BX , i);		\
      M3(flux, BX + Y , I,J,K, p) = F1(F_1d, BY , i);		\
      M3(flux, BX + Z , I,J,K, p) = F1(F_1d, BZ , i);		\
    }								\
} while (0)

  if (dir == 0) {
    PUT_LINE(0, 1, 2, i, j, k);
  } else if (dir == 1) {
    PUT_LINE(1, 2, 0, k, i, j);
  } else if (dir == 2) {
    PUT_LINE(2, 0, 1, j, k, i);
  } else {
    assert(0);
  }
#undef PUT_LINE
}

// ----------------------------------------------------------------------
// put_line_sc

// FIXME, make arg order consistent with put_line_fc
static void _mrc_unused
put_line_sc(struct mrc_fld *flux, struct mrc_fld *F,
	    int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PUT_LINE(X,Y,Z, I,J,K) do {					\
    for (int i = -l; i < ldim + r; i++) {				\
      M3(flux, RR   , I,J,K, p) = F1(F, RR , i);			\
      M3(flux, RVX+X, I,J,K, p) = F1(F, RVX, i);			\
      M3(flux, RVX+Y, I,J,K, p) = F1(F, RVY, i);			\
      M3(flux, RVX+Z, I,J,K, p) = F1(F, RVZ, i);			\
      M3(flux, UU   , I,J,K, p) = F1(F, UU , i);			\
    }									\
  } while (0)

  if (dir == 0) {
    PUT_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    PUT_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    PUT_LINE(2,0,1, j,k,i);
  }
#undef PUT_LINE
}

