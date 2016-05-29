
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

