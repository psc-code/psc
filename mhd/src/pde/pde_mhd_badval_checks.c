
#ifndef PDE_MHD_BADVAL_CHECKS_C
#define PDE_MHD_BADVAL_CHECKS_C

// FIXME, this is an extremely ugly way to use the same function for 
// the everything-in-one-field case and the separate-fields
// case

#ifndef GGCM_MHD_DEFS_EXTRA_H
#define _PP PP
#define _RR RR
#endif

// ----------------------------------------------------------------------
// patch_badval_checks_sc

static void
patch_badval_checks_sc(fld3d_t p_U, fld3d_t p_W)
{
  if (!s_do_badval_checks) {
    return;
  }

  int has_badval = 0;
  
  mrc_fld_data_t ppmin = 0.;
  mrc_fld_data_t rrmin = 0.;  // mhd->par.rrmin / mhd->rrnorm
  
  fld3d_foreach(i,j,k, 0, 0) {
    // check for negative pressure
    if (F3S(p_W, _PP, i,j,k) < ppmin) {
      has_badval = 5;
      mprintf("pressure @ (x=%g y=%g z=%g) = %lg < %lg\n",
	      PDE_CRDX_CC(i), PDE_CRDY_CC(j), PDE_CRDZ_CC(k), F3S(p_W, _PP, i,j,k), ppmin);
    }
    
    // check for negative density
    if (F3S(p_U, _RR, i,j,k) < rrmin) {
      has_badval = 4;
      mprintf("density @ (x=%g y=%g z=%g) = %lg < %lg\n",
	      PDE_CRDX_CC(i), PDE_CRDY_CC(j), PDE_CRDZ_CC(k), F3S(p_U, _RR, i,j,k), rrmin);
    }
    
    // check for invalid values
    for (int m = 0; m < 8; m++) {
      if (!isfinite(F3S(p_U, m, i,j,k))) {
	has_badval = 3;
	mprintf("NaN in field %d @ (x=%g y=%g z=%g)\n",
		m, PDE_CRDX_CC(i), PDE_CRDY_CC(j), PDE_CRDZ_CC(k));
      }
    }
  } fld3d_foreach_end;
  
  assert(!has_badval);
}

#endif
