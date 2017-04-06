
#ifndef PDE_MHD_BPUSH_C
#define PDE_MHD_BPUSH_C

// ----------------------------------------------------------------------
// patch_bpush1_c
//
// NOTE: this one expects -E, not E as input

static void
patch_bpush1_c(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev, fld3d_t p_E)
{
  if (p_Unext.arr_off == p_Uprev.arr_off) {
    fld3d_foreach(i,j,k, 0, 0) {
      F3S(p_Unext, BX, i,j,k) += dt * CURLX_FC(p_E, i,j,k);
      F3S(p_Unext, BY, i,j,k) += dt * CURLY_FC(p_E, i,j,k);
      F3S(p_Unext, BZ, i,j,k) += dt * CURLZ_FC(p_E, i,j,k);
    } fld3d_foreach_end;
  } else {
    fld3d_foreach(i,j,k, 0, 0) {
      F3S(p_Unext, BX, i,j,k) = F3S(p_Uprev, BX, i,j,k) + dt * CURLX_FC(p_E, i,j,k);
      F3S(p_Unext, BY, i,j,k) = F3S(p_Uprev, BY, i,j,k) + dt * CURLY_FC(p_E, i,j,k);
      F3S(p_Unext, BZ, i,j,k) = F3S(p_Uprev, BZ, i,j,k) + dt * CURLZ_FC(p_E, i,j,k);
    } fld3d_foreach_end;
  }
}

// ----------------------------------------------------------------------
// patch_bpush1_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define bpush1_F77 F77_FUNC(bpush1,BPUSH1)

void bpush1_F77(real *bx1, real *by1, real *bz1,
		real *bx3, real *by3, real *bz3,
		real *flx, real *fly, real *flz, real *dt);

static void
patch_bpush1_fortran(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev, fld3d_t p_E)
{
  bpush1_F77(F(p_Uprev, BX), F(p_Uprev, BY), F(p_Uprev, BZ), 
	     F(p_Unext, BX), F(p_Unext, BY), F(p_Unext, BZ), 
	     F(p_E, 0), F(p_E, 1), F(p_E, 2), &dt);
}

#endif

// ----------------------------------------------------------------------
// patch_bpush1

static void _mrc_unused
patch_bpush1(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev, fld3d_t p_E)
{
  if (s_opt_mhd_bpush1 == OPT_MHD_C) {
    patch_bpush1_c(p_Unext, dt, p_Uprev, p_E);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_bpush1 == OPT_MHD_FORTRAN) {
    patch_bpush1_fortran(p_Unext, dt, p_Uprev, p_E);
#endif
  } else {
    assert(0);
  }
}

#endif
