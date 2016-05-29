
// ----------------------------------------------------------------------
// pick_line_fc_cc

static void _mrc_unused
pick_line_fc_cc(fld1d_state_t u, struct mrc_fld *U,
		int ldim, int l, int r, int j, int k, int dim, int p)
{
#define PICK_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = -l; i < ldim + r; i++) {		\
      F1S(u, RR , i) = M3(U, RR   , I,J,K, p);		\
      F1S(u, RVX, i) = M3(U, RVX+X, I,J,K, p);		\
      F1S(u, RVY, i) = M3(U, RVX+Y, I,J,K, p);		\
      F1S(u, RVZ, i) = M3(U, RVX+Z, I,J,K, p);		\
      F1S(u, EE , i) = M3(U, EE   , I,J,K, p);		\
      F1S(u, BX , i) = M3(U, BX+X , I,J,K, p);		\
      F1S(u, BY , i) = M3(U, BX+Y , I,J,K, p);		\
      F1S(u, BZ , i) = M3(U, BX+Z , I,J,K, p);		\
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
// put_line_fc_cc

static void _mrc_unused
put_line_fc_cc(struct mrc_fld *flux, fld1d_state_t F,
	       int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = -l; i < ldim + r; i++) {			\
      M3(flux, RR     , I,J,K, p) = F1S(F, RR , i);		\
      M3(flux, RVX + X, I,J,K, p) = F1S(F, RVX, i);		\
      M3(flux, RVX + Y, I,J,K, p) = F1S(F, RVY, i);		\
      M3(flux, RVX + Z, I,J,K, p) = F1S(F, RVZ, i);		\
      M3(flux, EE     , I,J,K, p) = F1S(F, EE , i);		\
      M3(flux, BX + X , I,J,K, p) = F1S(F, BX , i);		\
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
// pick_line_fc

static void _mrc_unused
pick_line_fc(fld1d_state_t u
, struct mrc_fld *Bxi,
	     struct mrc_fld *U, struct mrc_fld *Bcc,
	     int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PICK_LINE(X, Y, Z, I, J, K)			       \
  do {							       \
    for (int i = - l; i < ldim + r; i++) {		       \
      F1S(u, RR , i) = M3(U, RR     , I,J,K, p);	       \
      F1S(u, RVX, i) = M3(U, RVX + X, I,J,K, p);	       \
      F1S(u, RVY, i) = M3(U, RVX + Y, I,J,K, p);	       \
      F1S(u, RVZ, i) = M3(U, RVX + Z, I,J,K, p);	       \
      F1S(u, EE , i) = M3(U, EE     , I,J,K, p);	       \
      F1S(u, BX , i) = M3(Bcc, X    , I,J,K, p);	       \
      F1S(u, BY , i) = M3(Bcc, Y    , I,J,K, p);	       \
      F1S(u, BZ , i) = M3(Bcc, Z    , I,J,K, p);	       \
      F1(Bxi, 0, i)  = M3(U, BX+X   , I,J,K, p);	       \
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
// put_line_fc

static void _mrc_unused
put_line_fc(struct mrc_fld *flux, fld1d_state_t F,
	    int ldim, int l, int r, int j, int k, int dir, int p)
{
#define PUT_LINE(X, Y, Z, I, J, K) do {				\
    for (int i = -l; i < ldim + r; i++) {			\
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
// pick_line_sc

static void _mrc_unused
pick_line_sc(fld1d_state_t u, struct mrc_fld *U,
	     int ldim, int l, int r, int j, int k, int dim, int p)
{
#define PICK_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = -l; i < ldim + r; i++) {		\
      F1S(u, RR , i) = M3(U, RR   , I,J,K, p);		\
      F1S(u, RVX, i) = M3(U, RVX+X, I,J,K, p);		\
      F1S(u, RVY, i) = M3(U, RVX+Y, I,J,K, p);		\
      F1S(u, RVZ, i) = M3(U, RVX+Z, I,J,K, p);		\
      F1S(u, UU , i) = M3(U, UU   , I,J,K, p);		\
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
// put_line_sc

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

