
#ifndef PDE_MHD_LINE_C
#define PDE_MHD_LINE_C

#include "pde/pde_mhd_setup.c"

// ----------------------------------------------------------------------
// mhd_get_line_state_fcons

static void _mrc_unused
mhd_get_line_state_fcons(fld1d_state_t u, struct mrc_fld *U,
			 int j, int k, int dir, int p, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1S(u, RR , i) = M3(U, RR   , I,J,K, p);		\
      F1S(u, RVX, i) = M3(U, RVX+X, I,J,K, p);		\
      F1S(u, RVY, i) = M3(U, RVX+Y, I,J,K, p);		\
      F1S(u, RVZ, i) = M3(U, RVX+Z, I,J,K, p);		\
      F1S(u, EE , i) = M3(U, EE   , I,J,K, p);		\
      F1S(u, BX , i) = M3(U, BX+X , I,J,K, p);		\
      F1S(u, BY , i) = M3(U, BX+Y , I,J,K, p);		\
      F1S(u, BZ , i) = M3(U, BX+Z , I,J,K, p);		\
      if (s_opt_divb == OPT_DIVB_GLM) {			\
	F1S(u, PSI, i) = M3(U, PSI, I,J,K, p);		\
      }							\
    }							\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_state_fcons

static void _mrc_unused
mhd_line_get_state_fcons(fld1d_state_t u, fld3d_t U,
			 int j, int k, int dir, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1S(u, RR , i) = F3S(U, RR   , I,J,K);		\
      F1S(u, RVX, i) = F3S(U, RVX+X, I,J,K);		\
      F1S(u, RVY, i) = F3S(U, RVX+Y, I,J,K);		\
      F1S(u, RVZ, i) = F3S(U, RVX+Z, I,J,K);		\
      F1S(u, EE , i) = F3S(U, EE   , I,J,K);		\
      F1S(u, BX , i) = F3S(U, BX+X , I,J,K);		\
      F1S(u, BY , i) = F3S(U, BX+Y , I,J,K);		\
      F1S(u, BZ , i) = F3S(U, BX+Z , I,J,K);		\
      if (s_opt_divb == OPT_DIVB_GLM) {			\
	F1S(u, PSI, i) = F3S(U, PSI, I,J,K);		\
      }							\
    }							\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_put_line_state_fcons

static void _mrc_unused
mhd_put_line_state_fcons(struct mrc_fld *flux, fld1d_state_t F,
			 int j, int k, int dir, int p, int ib, int ie)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = ib; i < ie; i++) {				\
      M3(flux, RR     , I,J,K, p) = F1S(F, RR , i);		\
      M3(flux, RVX + X, I,J,K, p) = F1S(F, RVX, i);		\
      M3(flux, RVX + Y, I,J,K, p) = F1S(F, RVY, i);		\
      M3(flux, RVX + Z, I,J,K, p) = F1S(F, RVZ, i);		\
      M3(flux, EE     , I,J,K, p) = F1S(F, EE , i);		\
      M3(flux, BX + X , I,J,K, p) = F1S(F, BX , i);		\
      M3(flux, BX + Y , I,J,K, p) = F1S(F, BY , i);		\
      M3(flux, BX + Z , I,J,K, p) = F1S(F, BZ , i);		\
      if (s_opt_divb == OPT_DIVB_GLM) {				\
	M3(flux, PSI, I,J,K, p)   = F1S(F, PSI, i);		\
      }								\
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
// mhd_line_put_state_fcons

static void _mrc_unused
mhd_line_put_state_fcons(fld1d_state_t l_U, fld3d_t p_U,
			 int j, int k, int dir, int ib, int ie)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {					\
    for (int i = ib; i < ie; i++) {					\
      F3S(p_U, RR     , I,J,K) = F1S(l_U, RR , i);			\
      F3S(p_U, RVX + X, I,J,K) = F1S(l_U, RVX, i);			\
      F3S(p_U, RVX + Y, I,J,K) = F1S(l_U, RVY, i);			\
      F3S(p_U, RVX + Z, I,J,K) = F1S(l_U, RVZ, i);			\
      F3S(p_U, EE     , I,J,K) = F1S(l_U, EE , i);			\
      F3S(p_U, BX + X , I,J,K) = F1S(l_U, BX , i);			\
      F3S(p_U, BX + Y , I,J,K) = F1S(l_U, BY , i);			\
      F3S(p_U, BX + Z , I,J,K) = F1S(l_U, BZ , i);			\
      if (s_opt_divb == OPT_DIVB_GLM) {					\
	F3S(p_U, PSI, I,J,K)   = F1S(l_U, PSI, i);			\
      }									\
    }									\
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
// mhd_get_line_state_fcons_ct

static void _mrc_unused
mhd_get_line_state_fcons_ct(fld1d_state_t u, fld1d_t bx,
			    struct mrc_fld *U, struct mrc_fld *Bcc,
			    int j, int k, int dir, int p, int ib, int ie)
{
#define GET_LINE(X, Y, Z, I, J, K)			       \
  do {							       \
    for (int i = ib; i < ie; i++) {			       \
      F1S(u, RR , i) = M3(U, RR     , I,J,K, p);	       \
      F1S(u, RVX, i) = M3(U, RVX + X, I,J,K, p);	       \
      F1S(u, RVY, i) = M3(U, RVX + Y, I,J,K, p);	       \
      F1S(u, RVZ, i) = M3(U, RVX + Z, I,J,K, p);	       \
      F1S(u, EE , i) = M3(U, EE     , I,J,K, p);	       \
      F1S(u, BX , i) = M3(Bcc, X    , I,J,K, p);	       \
      F1S(u, BY , i) = M3(Bcc, Y    , I,J,K, p);	       \
      F1S(u, BZ , i) = M3(Bcc, Z    , I,J,K, p);	       \
      F1(bx, i)      = M3(U, BX+X   , I,J,K, p);	       \
    }							       \
  } while (0)
  if (dir == 0) {
    GET_LINE(0, 1, 2, i, j, k);
  } else if (dir == 1) {
    GET_LINE(1, 2, 0, k, i, j);
  } else if (dir == 2) {
    GET_LINE(2, 0, 1, j, k, i);
  } else {
    assert(0);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_put_line_state_fcons_ct

static void _mrc_unused
mhd_put_line_state_fcons_ct(struct mrc_fld *flux, fld1d_state_t F,
			    int j, int k, int dir, int p, int ib, int ie)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = ib; i < ie; i++) {				\
      M3(flux, RR     , I,J,K, p) = F1S(F, RR , i);		\
      M3(flux, RVX + X, I,J,K, p) = F1S(F, RVX, i);		\
      M3(flux, RVX + Y, I,J,K, p) = F1S(F, RVY, i);		\
      M3(flux, RVX + Z, I,J,K, p) = F1S(F, RVZ, i);		\
      M3(flux, EE     , I,J,K, p) = F1S(F, EE , i);		\
      M3(flux, BX + Y , I,J,K, p) = F1S(F, BY , i);		\
      M3(flux, BX + Z , I,J,K, p) = F1S(F, BZ , i);		\
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
// mhd_get_line_state_scons

static void _mrc_unused
mhd_get_line_state_scons(fld1d_state_t u, struct mrc_fld *U,
			 int j, int k, int dim, int p, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1S(u, RR , i) = M3(U, RR   , I,J,K, p);		\
      F1S(u, RVX, i) = M3(U, RVX+X, I,J,K, p);		\
      F1S(u, RVY, i) = M3(U, RVX+Y, I,J,K, p);		\
      F1S(u, RVZ, i) = M3(U, RVX+Z, I,J,K, p);		\
      F1S(u, UU , i) = M3(U, UU   , I,J,K, p);		\
    }							\
  } while (0)

  if (dim == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dim == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dim == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_state_scons

static void _mrc_unused
mhd_line_get_state_scons(fld1d_state_t u, fld3d_t U,
			 int j, int k, int dim, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1S(u, RR , i) = F3S(U, RR   , I,J,K);		\
      F1S(u, RVX, i) = F3S(U, RVX+X, I,J,K);		\
      F1S(u, RVY, i) = F3S(U, RVX+Y, I,J,K);		\
      F1S(u, RVZ, i) = F3S(U, RVX+Z, I,J,K);		\
      F1S(u, UU , i) = F3S(U, UU   , I,J,K);		\
    }							\
  } while (0)

  if (dim == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dim == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dim == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_put_line_state_scons

static void _mrc_unused
mhd_put_line_state_scons(struct mrc_fld *flux, fld1d_state_t F,
			 int j, int k, int dir, int p, int ib, int ie)
{
#define PUT_LINE(X,Y,Z, I,J,K) do {					\
    for (int i = ib; i < ie; i++) {					\
      M3(flux, RR   , I,J,K, p) = F1S(F, RR , i);			\
      M3(flux, RVX+X, I,J,K, p) = F1S(F, RVX, i);			\
      M3(flux, RVX+Y, I,J,K, p) = F1S(F, RVY, i);			\
      M3(flux, RVX+Z, I,J,K, p) = F1S(F, RVZ, i);			\
      M3(flux, UU   , I,J,K, p) = F1S(F, UU , i);			\
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

// ----------------------------------------------------------------------
// mhd_line_put_state_scons

static void _mrc_unused
mhd_line_put_state_scons(fld1d_state_t l_U, fld3d_t p_U,
			 int j, int k, int dir, int ib, int ie)
{
#define PUT_LINE(X,Y,Z, I,J,K) do {					\
    for (int i = ib; i < ie; i++) {					\
      F3S(p_U, RR   , I,J,K) = F1S(l_U, RR , i);			\
      F3S(p_U, RVX+X, I,J,K) = F1S(l_U, RVX, i);			\
      F3S(p_U, RVX+Y, I,J,K) = F1S(l_U, RVY, i);			\
      F3S(p_U, RVX+Z, I,J,K) = F1S(l_U, RVZ, i);			\
      F3S(p_U, UU   , I,J,K) = F1S(l_U, UU , i);			\
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

// ----------------------------------------------------------------------
// mhd_get_line_state

static void _mrc_unused
mhd_get_line_state(fld1d_state_t u, struct mrc_fld *U, int j, int k, int dir,
                  int p, int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_get_line_state_fcons(u, U, j, k, dir, p, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_get_line_state_scons(u, U, j, k, dir, p, ib, ie);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_line_get_state

static void _mrc_unused
mhd_line_get_state(fld1d_state_t u, fld3d_t U, int j, int k, int dir,
		   int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_line_get_state_fcons(u, U, j, k, dir, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_line_get_state_scons(u, U, j, k, dir, ib, ie);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_put_line_state

static void _mrc_unused
mhd_put_line_state(struct mrc_fld *U, fld1d_state_t u,
		   int j, int k, int dir, int p, int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_put_line_state_fcons(U, u, j, k, dir, p, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_put_line_state_scons(U, u, j, k, dir, p, ib, ie);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_line_put_state

static void _mrc_unused
mhd_line_put_state(fld1d_state_t l_U, fld3d_t p_U,
		   int j, int k, int dir, int ib, int ie)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mhd_line_put_state_fcons(l_U, p_U, j, k, dir, ib, ie);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mhd_line_put_state_scons(l_U, p_U, j, k, dir, ib, ie);
  } else {
    assert(0);
  }
}


// ----------------------------------------------------------------------
// mhd_get_line_1

static void _mrc_unused
mhd_get_line_1(fld1d_t u, struct mrc_fld *U,
	       int j, int k, int dir, int p, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1(u, i) = M3(U, 0, I,J,K, p);			\
    }							\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_1

static void _mrc_unused
mhd_line_get_1(fld1d_t u, fld3d_t U,
	       int j, int k, int dir, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1(u, i) = F3S(U, 0, I,J,K);			\
    }							\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_vec

static void _mrc_unused
mhd_line_get_vec(fld1d_vec_t u, fld3d_t U,
		 int j, int k, int dir, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = ib; i < ie; i++) {			\
      F1S(u, 0, i) = F3S(U, X , I,J,K);			\
      F1S(u, 1, i) = F3S(U, Y , I,J,K);			\
      F1S(u, 2, i) = F3S(U, Z , I,J,K);			\
    }							\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_vec_fc
// 
// get 1d line of values from c.c. fld3d_t, average to faces

static void _mrc_unused
mhd_line_get_vec_fc(fld1d_vec_t u, fld3d_t U,
		    int j, int k, int dir, int ib, int ie)
{
#define GET_LINE(X,Y,Z,I,J,K,DI,DJ,DK) do {				\
    for (int i = ib; i < ie; i++) {					\
      F1S(u, 0, i) = .5f*(F3S(U, X , I,J,K) + F3S(U, X, I-DI,J-DJ,K-DK)); \
      F1S(u, 1, i) = .5f*(F3S(U, Y , I,J,K) + F3S(U, Y, I-DI,J-DJ,K-DK)); \
      F1S(u, 2, i) = .5f*(F3S(U, Z , I,J,K) + F3S(U, Z, I-DI,J-DJ,K-DK)); \
    }									\
  } while (0)

  if (dir == 0) {
    GET_LINE(0,1,2, i,j,k, 1,0,0);
  } else if (dir == 1) {
    GET_LINE(1,2,0, k,i,j, 0,1,0);
  } else if (dir == 2) {
    GET_LINE(2,0,1, j,k,i, 0,0,1);
  }
#undef GET_LINE
}

// ----------------------------------------------------------------------
// mhd_line_get_b0

static void _mrc_unused
mhd_line_get_b0(int j, int k, int dir, int ib, int ie)
{
  if (!s_opt_background) {
    return;
  }

  mhd_line_get_vec(s_aux.b0, s_p_aux.b0, j, k, dir, ib, ie + 1);
}

// ----------------------------------------------------------------------
// mhd_line_get_b0_fc

static void _mrc_unused
mhd_line_get_b0_fc(int j, int k, int dir, int ib, int ie)
{
  if (!s_opt_background) {
    return;
  }

  mhd_line_get_vec_fc(s_aux.b0, s_p_aux.b0, j, k, dir, ib, ie + 1);
}

#endif
