
// ----------------------------------------------------------------------
// face-centered current macros
//
// D0/1/2 refers to which face

// FIXME, it'd be nice to have a mechanism to not calculate the zero terms in < 3 dims

// FIXME, s_p_inv_dx* needs shifting, too

#define D0_x(x, m, i,j,k, p) ((M3(x, m, i   ,j,k, p) - M3(x, m, i-di,j,k, p)) * PDE_INV_DX(i))
#define D0_y(x, m, i,j,k, p) (((M3(x, m, i-di,j+dj,k, p) + M3(x, m, i,j+dj,k, p)) - \
			       (M3(x, m, i-di,j-dj,k, p) + M3(x, m, i,j-dj,k, p))) * .25f * PDE_INV_DYF(j))
#define D0_z(x, m, i,j,k, p) (((M3(x, m, i-di,j,k+dk, p) + M3(x, m, i,j,k+dk, p)) - \
			       (M3(x, m, i-di,j,k-dk, p) + M3(x, m, i,j,k-dk, p))) * .25f * PDE_INV_DZF(k))

#define D1_x(x, m, i,j,k, p) (((M3(x, m, i+di,j-dj,k, p) + M3(x, m, i+di,j,k, p)) - \
			       (M3(x, m, i-di,j-dj,k, p) + M3(x, m, i-di,j,k, p))) * .25f * PDE_INV_DXF(i))
#define D1_y(x, m, i,j,k, p) ((M3(x, m, i,j   ,k, p) - M3(x, m, i,j-dj,k, p)) * PDE_INV_DY(j))
#define D1_z(x, m, i,j,k, p) (((M3(x, m, i,j-dj,k+dk, p) + M3(x, m, i,j,k+dk, p)) - \
			       (M3(x, m, i,j-dj,k-dk, p) + M3(x, m, i,j,k-dk, p))) * .25f * PDE_INV_DZF(k))

#define D2_x(x, m, i,j,k, p) (((M3(x, m, i+di,j,k-dk, p) + M3(x, m, i+di,j,k, p)) - \
			       (M3(x, m, i-di,j,k-dk, p) + M3(x, m, i-di,j,k, p))) * .25f * PDE_INV_DXF(i))
#define D2_y(x, m, i,j,k, p) (((M3(x, m, i,j+dj,k-dk, p) + M3(x, m, i,j+dj,k, p)) - \
			       (M3(x, m, i,j-dj,k-dk, p) + M3(x, m, i,j-dj,k, p))) * .25f * PDE_INV_DYF(j))
#define D2_z(x, m, i,j,k, p) ((M3(x, m, i,j,k   , p) - M3(x, m, i,j,k-dk, p)) * PDE_INV_DZ(k))

// ----------------------------------------------------------------------
// calc_J_line

static void _mrc_unused
calc_J_line(fld1d_vec_t J, struct mrc_fld *x,
	    int jj, int kk, int dir, int p, int ib, int ie)
{
  if (dir == 0) {
    int j = jj, k = kk;
    for (int i = ib; i < ie; i++) {
      F1V(J, 0, i) = D0_y(x, BZ, i,j,k, p) - D0_z(x, BY, i,j,k, p);
      F1V(J, 1, i) = D0_z(x, BX, i,j,k, p) - D0_x(x, BZ, i,j,k, p);
      F1V(J, 2, i) = D0_x(x, BY, i,j,k, p) - D0_y(x, BX, i,j,k, p);
    }
  } else if (dir == 1) {
    int k = jj, i = kk;
    for (int j = ib; j < ie; j++) {
      F1V(J, 2, j) = D1_y(x, BZ, i,j,k, p) - D1_z(x, BY, i,j,k, p);
      F1V(J, 0, j) = D1_z(x, BX, i,j,k, p) - D1_x(x, BZ, i,j,k, p);
      F1V(J, 1, j) = D1_x(x, BY, i,j,k, p) - D1_y(x, BX, i,j,k, p);
    }
  } else if (dir == 2) {
    int i = jj, j = kk;
    for (int k = ib; k < ie; k++) {
      F1V(J, 1, k) = D2_y(x, BZ, i,j,k, p) - D2_z(x, BY, i,j,k, p);
      F1V(J, 2, k) = D2_z(x, BX, i,j,k, p) - D2_x(x, BZ, i,j,k, p);
      F1V(J, 0, k) = D2_x(x, BY, i,j,k, p) - D2_y(x, BX, i,j,k, p);
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_line_get_current

static void _mrc_unused
mhd_line_get_current(struct mrc_fld *x, int j, int k, int dir,
		     int p, int ib, int ie)
{
  if (!s_opt_need_current) {
    return;
  }

  calc_J_line(s_aux.j, x, j, k, dir, p, ib, ie);
}

