
// ----------------------------------------------------------------------
// pick_line_fc_cc

static void _mrc_unused
pick_line_fc_cc(struct mrc_fld *x1, struct mrc_fld *x,
		int ldim, int l, int r, int j, int k, int dim, int p)
{
#define PICK_LINE(X,Y,Z,I,J,K) do {			\
    for (int i = -l; i < ldim + r; i++) {		\
      F1(x1, RR , i) = M3(x, RR   , I,J,K, p);		\
      F1(x1, RVX, i) = M3(x, RVX+X, I,J,K, p);		\
      F1(x1, RVY, i) = M3(x, RVX+Y, I,J,K, p);		\
      F1(x1, RVZ, i) = M3(x, RVX+Z, I,J,K, p);		\
      F1(x1, EE , i) = M3(x, EE   , I,J,K, p);		\
      F1(x1, BX , i) = M3(x, BX+X , I,J,K, p);		\
      F1(x1, BY , i) = M3(x, BX+Y , I,J,K, p);		\
      F1(x1, BZ , i) = M3(x, BX+Z , I,J,K, p);		\
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

