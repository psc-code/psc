
// ----------------------------------------------------------------------
// face-centered current macros
//
// D0/1/2 refers to which face

// FIXME, it'd be nice to have a mechanism to not calculate the zero terms in < 3 dims

// FIXME, s_p_inv_dx* needs shifting, too

#define D0_x(x, m, i,j,k) ((F3S(x, m, i   ,j,k) - F3S(x, m, i-di,j,k)) * PDE_INV_DX(i))
#define D0_y(x, m, i,j,k) (((F3S(x, m, i-di,j+dj,k) + F3S(x, m, i,j+dj,k)) - \
			    (F3S(x, m, i-di,j-dj,k) + F3S(x, m, i,j-dj,k))) * .25f * PDE_INV_DYF(j))
#define D0_z(x, m, i,j,k) (((F3S(x, m, i-di,j,k+dk) + F3S(x, m, i,j,k+dk)) - \
			    (F3S(x, m, i-di,j,k-dk) + F3S(x, m, i,j,k-dk))) * .25f * PDE_INV_DZF(k))

#define D1_x(x, m, i,j,k) (((F3S(x, m, i+di,j-dj,k) + F3S(x, m, i+di,j,k)) - \
			    (F3S(x, m, i-di,j-dj,k) + F3S(x, m, i-di,j,k))) * .25f * PDE_INV_DXF(i))
#define D1_y(x, m, i,j,k) ((F3S(x, m, i,j   ,k) - F3S(x, m, i,j-dj,k)) * PDE_INV_DY(j))
#define D1_z(x, m, i,j,k) (((F3S(x, m, i,j-dj,k+dk) + F3S(x, m, i,j,k+dk)) - \
			    (F3S(x, m, i,j-dj,k-dk) + F3S(x, m, i,j,k-dk))) * .25f * PDE_INV_DZF(k))

#define D2_x(x, m, i,j,k) (((F3S(x, m, i+di,j,k-dk) + F3S(x, m, i+di,j,k)) - \
			    (F3S(x, m, i-di,j,k-dk) + F3S(x, m, i-di,j,k))) * .25f * PDE_INV_DXF(i))
#define D2_y(x, m, i,j,k) (((F3S(x, m, i,j+dj,k-dk) + F3S(x, m, i,j+dj,k)) - \
			    (F3S(x, m, i,j-dj,k-dk) + F3S(x, m, i,j-dj,k))) * .25f * PDE_INV_DYF(j))
#define D2_z(x, m, i,j,k) ((F3S(x, m, i,j,k   ) - F3S(x, m, i,j,k-dk)) * PDE_INV_DZ(k))

// ----------------------------------------------------------------------
// calc_J_line

static void _mrc_unused
calc_J_line(fld1d_vec_t J, fld3d_t x,
	    int jj, int kk, int dir, int ib, int ie)
{
  if (dir == 0) {
    int j = jj, k = kk;
    for (int i = ib; i < ie; i++) {
      F1V(J, 0, i) = D0_y(x, BZ, i,j,k) - D0_z(x, BY, i,j,k);
      F1V(J, 1, i) = D0_z(x, BX, i,j,k) - D0_x(x, BZ, i,j,k);
      F1V(J, 2, i) = D0_x(x, BY, i,j,k) - D0_y(x, BX, i,j,k);
    }
  } else if (dir == 1) {
    int k = jj, i = kk;
    for (int j = ib; j < ie; j++) {
      F1V(J, 2, j) = D1_y(x, BZ, i,j,k) - D1_z(x, BY, i,j,k);
      F1V(J, 0, j) = D1_z(x, BX, i,j,k) - D1_x(x, BZ, i,j,k);
      F1V(J, 1, j) = D1_x(x, BY, i,j,k) - D1_y(x, BX, i,j,k);
    }
  } else if (dir == 2) {
    int i = jj, j = kk;
    for (int k = ib; k < ie; k++) {
      F1V(J, 1, k) = D2_y(x, BZ, i,j,k) - D2_z(x, BY, i,j,k);
      F1V(J, 2, k) = D2_z(x, BX, i,j,k) - D2_x(x, BZ, i,j,k);
      F1V(J, 0, k) = D2_x(x, BY, i,j,k) - D2_y(x, BX, i,j,k);
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_line_get_current

static void _mrc_unused
mhd_line_get_current(fld3d_t x, int j, int k, int dir,
		     int ib, int ie)
{
  if (!s_opt_need_current) {
    return;
  }

  calc_J_line(s_aux.j, x, j, k, dir, ib, ie);
}

