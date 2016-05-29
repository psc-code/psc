
#include "ggcm_mhd_private.h"

#include <math.h>

#define TINY_NUMBER 1.0e-20 // FIXME

// ======================================================================
// semi-conservative

// ======================================================================

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

