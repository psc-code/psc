
#define BOUNDS_CHECK

#include "ggcm_mhd_bndsw.h"
#include "ggcm_mhd_crds.h"
#include "mrc_domain.h"
#include "ggcm_mhd_gkeyll.h"

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

struct ggcm_mhd_bnd_sub {
  double bnvals[N_PRIMITIVE];
  bool apply_bndsw;
  bool do_legacy;
};

#define ggcm_mhd_bnd_sub(bnd) mrc_to_subobj(bnd, struct ggcm_mhd_bnd_sub);

// FIXME
// The B boundary conditions have numerous issues:
//
// - It's not clear how many ghost points we need in the first place, and how
//   this all interacts with the bpush update (I think bpush updates a boundary
//   point, which however afterwards will be fixed up here, so should be okay)
//   The current implementation of eta J in Ohm's Law probably really needs two
//   ghost points.
//   There is still asymmetric behavior, e.g., setting normal B field in obndra_zh
//   is required to avoid things going wrong at the zh-xh edge, while there's no
//   problem at the zl-xh edge. working guess is that this is related to the pusher,
//   actually
//
// - The "Set div B = 0" implementation at corners is iffy (order dependent), it
//   will, e.g., set By based on Bx at the wall (which is not well defined), and
//   then Bx at the wall based on that By. It's unlikely to really matter, though.
//   (But it might be the reason for the remaining small discrepancies between the
//   two staggerings.)
//
// - Known difference to the Fortran version: Bnormal components in obndra() are copied
//   with their staggering in mind now, rather than treating them as if they were c.c.
//
// - We really should only be using _BX/_BY/_BZ, but I think they're used in all cases
//   that actually matter -- changing the other places would need
//   adjusting of the loop limits, while now we just set one more, unneeded point on either side
//   depending on sc_ggcm / sc (But that's just a guess)

// FIXME, which of these ghost points are actually used? / loop limits

#define _BX(f, ix,iy,iz, p) M3(f, BX, ix+SHIFT,iy,iz, p)
#define _BY(f, ix,iy,iz, p) M3(f, BY, ix,iy+SHIFT,iz, p)
#define _BZ(f, ix,iy,iz, p) M3(f, BZ, ix,iy,iz+SHIFT, p)

#if MT == MT_FCONS_CC
// maybe this should use fd1 to be more explicit, but that's identical to bd3, anyway

#define BNDDIV_BY_L(ix, iy, iz, p)					\
  (_BY(f, ix,iy+2,iz, p) +						\
   bdx3[ix]/bdy3[iy+1] * (_BX(f, ix+1,iy+1,iz  , p) -			\
			  _BX(f, ix-1,iy+1,iz  , p)) +			\
   bdz3[iz]/bdy3[iy+1] * (_BZ(f, ix  ,iy+1,iz+1, p) -			\
			  _BZ(f, ix  ,iy+1,iz-1, p)))
#define BNDDIV_BZ_L(ix, iy, iz, p)					\
  (_BZ(f, ix,iy,iz+2, p) +						\
   bdx3[ix]/bdz3[iz+1] * (_BX(f, ix+1,iy  ,iz+1, p) -			\
			  _BX(f, ix-1,iy  ,iz+1, p)) +			\
   bdy3[iy]/bdz3[iz+1] * (_BY(f, ix  ,iy+1,iz+1, p) -			\
			  _BY(f, ix  ,iy-1,iz+1, p)))

#define BNDDIV_BX_H(ix, iy, iz, p)					\
  (_BX(f, ix-2,iy,iz, p) -						\
   bdy3[iy]/bdx3[ix-1] * (_BY(f, ix-1,iy+1,iz  , p) -			\
			  _BY(f, ix-1,iy-1,iz  , p)) -			\
   bdz3[iz]/bdx3[ix-1] * (_BZ(f, ix-1,iy  ,iz+1, p) -			\
			  _BZ(f, ix-1,iy  ,iz-1, p)))
#define BNDDIV_BY_H(ix, iy, iz, p)					\
  (_BY(f, ix,iy-2,iz, p) -						\
   bdx3[ix]/bdy3[iy-1] * (_BX(f, ix+1,iy-1,iz  , p) -			\
			  _BX(f, ix-1,iy-1,iz  , p)) -			\
   bdz3[iz]/bdy3[iy-1] * (_BZ(f, ix  ,iy-1,iz+1, p) -			\
			  _BZ(f, ix  ,iy-1,iz-1, p)))
#define BNDDIV_BZ_H(ix, iy, iz, p)					\
  (_BZ(f, ix,iy,iz-2, p) -						\
   bdx3[ix]/bdz3[iz-1] * (_BX(f, ix+1,iy  ,iz-1, p) -			\
			  _BX(f, ix-1,iy  ,iz-1, p)) -			\
   bdy3[iy]/bdz3[iz-1] * (_BY(f, ix  ,iy+1,iz-1, p) -			\
			  _BY(f, ix  ,iy-1,iz-1, p)))

#else

#define BNDDIV_BY_L(ix, iy, iz, p)					\
  (_BY(f, ix,iy+1,iz, p) +						\
   bdx3[ix]/bdy3[iy] * (_BX(f, ix+1,iy,iz  , p) -			\
			_BX(f, ix  ,iy,iz  , p)) +			\
   bdz3[iz]/bdy3[iy] * (_BZ(f, ix  ,iy,iz+1, p) -			\
			_BZ(f, ix  ,iy,iz  , p)))
#define BNDDIV_BZ_L(ix, iy, iz, p)					\
  (_BZ(f, ix,iy,iz+1, p) +						\
   bdx3[ix]/bdz3[iz] * (_BX(f, ix+1,iy  ,iz, p) -			\
			_BX(f, ix  ,iy  ,iz, p)) +			\
   bdy3[iy]/bdz3[iz] * (_BY(f, ix  ,iy+1,iz, p) -			\
			_BY(f, ix  ,iy  ,iz, p)))

#define BNDDIV_BX_H(ix, iy, iz, p)					\
  (_BX(f, ix-1,iy,iz, p) -						\
   bdy3[iy]/bdx3[ix-1] * (_BY(f, ix-1,iy+1,iz  , p) -			\
			  _BY(f, ix-1,iy  ,iz  , p)) -			\
   bdz3[iz]/bdx3[ix-1] * (_BZ(f, ix-1,iy  ,iz+1, p) -			\
			  _BZ(f, ix-1,iy  ,iz  , p)))
#define BNDDIV_BY_H(ix, iy, iz, p)					\
  (_BY(f, ix,iy-1,iz, p) -						\
   bdx3[ix]/bdy3[iy-1] * (_BX(f, ix+1,iy-1,iz  , p) -			\
			  _BX(f, ix  ,iy-1,iz  , p)) -			\
   bdz3[iz]/bdy3[iy-1] * (_BZ(f, ix  ,iy-1,iz+1, p) -			\
			  _BZ(f, ix  ,iy-1,iz  , p)))
#define BNDDIV_BZ_H(ix, iy, iz, p)					\
  (_BZ(f, ix,iy,iz-1, p) -						\
   bdx3[ix]/bdz3[iz-1] * (_BX(f, ix+1,iy  ,iz-1, p) -			\
			  _BX(f, ix  ,iy  ,iz-1, p)) -			\
   bdy3[iy]/bdz3[iz-1] * (_BY(f, ix  ,iy+1,iz-1, p) -			\
			  _BY(f, ix  ,iy  ,iz-1, p)))

#endif

// FIXME, mv -> mrc_domain
static bool
mrc_domain_at_boundary_lo(struct mrc_domain *domain, int d, int p)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(domain, p, &info);
  return (info.off[d] == 0);
}

static bool
mrc_domain_at_boundary_hi(struct mrc_domain *domain, int d, int p)
{
  int gdims[3]; mrc_domain_get_global_dims(domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(domain, p, &info);
  return (info.off[d] + info.ldims[d] == gdims[d] * (1 << info.level));
}

// ----------------------------------------------------------------------
// bnd_sw

static void
bnd_sw(struct ggcm_mhd_bnd *bnd, int ix, int iy, int iz, int p, float bn[], float bntim)
{
  struct ggcm_mhd_bnd_sub *sub = ggcm_mhd_bnd_sub(bnd);
  struct ggcm_mhd *mhd = bnd->mhd;

  static bool first_time = true;
  static struct ggcm_mhd_bndsw *bndsw;
  if (first_time) {
    bndsw = ggcm_mhd_get_var_obj(mhd, "bndsw");
  }
  first_time = false;

  if (bndsw) {
    struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
    float xx[3] = { MRC_MCRDX(crds, ix, p),
		    MRC_MCRDY(crds, iy, p),
		    MRC_MCRDZ(crds, iz, p), };

    ggcm_mhd_bndsw_at(bndsw, bntim, xx, bn);
  } else {
    bn[RR] = sub->bnvals[RR] / mhd->rrnorm;
    bn[VX] = sub->bnvals[VX] / mhd->vvnorm;
    bn[VY] = sub->bnvals[VY] / mhd->vvnorm;
    bn[VZ] = sub->bnvals[VZ] / mhd->vvnorm;
    bn[PP] = sub->bnvals[PP] / mhd->ppnorm;
    bn[BX] = sub->bnvals[BX] / mhd->bbnorm;
    bn[BY] = sub->bnvals[BY] / mhd->bbnorm;
    bn[BZ] = sub->bnvals[BZ] / mhd->bbnorm;
  }
}

// ----------------------------------------------------------------------
// obndra_xl_bndsw
//
// set inflow fluid boundary conditions at x-low boundary

static void
obndra_xl_bndsw(struct ggcm_mhd_bnd *bnd, struct mrc_fld *f, float bntim, int p)
{
  struct ggcm_mhd *mhd = bnd->mhd;
  struct mrc_fld *b0 = NULL;
  if (mhd->b0) {
    b0 = mrc_fld_get_as(mhd->b0, FLD_TYPE); // FIXME needed?
  }
  const int *sw = mrc_fld_spatial_sw(f), *dims = mrc_fld_spatial_dims(f);

  int swx = sw[0], swy = sw[1], swz = sw[2];
  int my = dims[1], mz = dims[2];

  for (int iz = -swz; iz < mz + swz; iz++) {
    for (int iy = -swy; iy < my + swy; iy++) {
      for (int ix = -swx; ix < 0; ix++) {
        float bn[N_PRIMITIVE];
        bnd_sw(bnd, ix, iy, iz, p, bn, bntim);

	// FIXME, the background field handling is no good.
	// There are two competing objectives: convert_state_from_prim() needs to know the B_1
	// in the fully conservative case to calculate internal energy correctly. OTOH, it needs to
	// know the full B = B_0 + B_1 to calculate the E = - v x B correctly in the gkeyll case.
	// Both cases are covered below, but obviously that's not pretty.
	// I think the solution is to pass B and B_0, and have the conversion do the right thing.
	
#if MT_FORMULATION(MT) != MT_FORMULATION_GKEYLL
	// subtract background field if used
	if (b0) {
	  for (int d = 0; d < 3; d++) {
	    bn[BX + d] -= M3(b0, d, ix,iy,iz, p);
	  }
	}
#endif

	mrc_fld_data_t prim[N_PRIMITIVE], state[s_n_state];
	for (int m = 0; m < N_PRIMITIVE; m++) {
	  prim[m] = bn[m];
	}
	convert_state_from_prim(state, prim);

	if (MT_BGRID(MT) == MT_BGRID_CC) {
	  convert_put_state_to_3d(state, f, ix,iy,iz, p);

#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
	  if (b0) {
	    for (int d = 0; d < 3; d++) {
	      M3(f, s_gk_idx_em + GK_BX + d, ix,iy,iz, p) -= M3(b0, d, ix,iy,iz, p);
	    }
	  }
#endif
	} else { // staggered B
	  convert_put_fluid_state_to_3d(state, f, ix,iy,iz, p);
	  _BX(f, ix+1,iy,iz, p) = state[BX];
	  if (iy > -swy) {
	    _BY(f, ix,iy,iz, p) = state[BY];
	  }
	  if (iz > -swz) {
	    _BZ(f, ix,iy,iz, p) = state[BZ];
	  }
	}
      }
    }
  }

  if (mhd->b0) {
    mrc_fld_put_as(b0, mhd->b0);
  }
}


// ----------------------------------------------------------------------
// obndra_yl_open

static void
obndra_yl_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       const int sw[3], const int ldims[3], int b, int p)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 1, p)) {
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iy = 0; iy < sw[1]; iy++) {
	  for (int m = 0; m < 5; m++) {
	    M3(f, m, ix,b-iy,iz, p) = M3(f, m, ix,b+1+iy,iz, p);
	  }
	  if (ix > -sw[0]) {
	    _BX(f, ix,b-iy,iz, p) = _BX(f, ix,b+1+iy,iz, p);
	  }
	  if (iz > -sw[2]) {
	    _BZ(f, ix,b-iy,iz, p) = _BZ(f, ix,b+1+iy,iz, p);
	  }
	}
	for (int iy = 1; iy < sw[1]; iy++) {
	  _BY(f, ix,-iy,iz, p) = _BY(f, ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_zl_open

static void
obndra_zl_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       const int sw[3], const int ldims[3], int b, int p)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 2, p)) {
    for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iz = 0; iz < sw[2]; iz++) {
	  for (int m = 0; m < 5; m++) {
	    M3(f, m, ix,iy,b-iz, p) = M3(f, m, ix,iy,b+1+iz, p);
	  }
	  if (ix > -sw[0]) {
	    _BX(f, ix,iy,b-iz, p) = _BX(f, ix,iy,b+1+iz, p);
	  }
	  if (iy > -sw[1]) {
	    _BY(f, ix,iy,b-iz, p) = _BY(f, ix,iy,b+1+iz, p);
	  }
	}
	for (int iz = 1; iz < sw[2]; iz++) {
	  _BZ(f, ix,iy,-iz, p) = _BZ(f, ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_xh_open

static void
obndra_xh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       const int sw[3], const int ldims[3], int b, int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 0, p)) {
    int mx = ldims[0];
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
	for (int ix = 0; ix < sw[0]; ix++) {
	  for (int m = 0; m < 5; m++) {
	    M3(f, m, b+ix,iy,iz, p) = M3(f, m, b-ix-1,iy,iz, p);
	  }
	  if (iy > -sw[1]) {
	    _BY(f, b+ix,iy,iz, p) = _BY(f, b-ix-1,iy,iz, p);
	  }
	  if (iz > -sw[2]) {
	    _BZ(f, b+ix,iy,iz, p) = _BZ(f, b-ix-1,iy,iz, p);
	  }
	}
	for (int ix = 1; ix < sw[0]; ix++) {
	  _BX(f, mx+ix,iy,iz, p) = _BX(f, mx-ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_yh_open

static void
obndra_yh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       const int sw[3], const int ldims[3], int b, int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 1, p)) {
    int my = ldims[1];
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iy = 0; iy < sw[1]; iy++) {
	  for (int m = 0; m < 5; m++) {
	    M3(f, m, ix,b+iy,iz, p) = M3(f,m, ix,b-iy-1,iz, p);
	  }
	  if (ix > -sw[0]) {
	    _BX(f, ix,b+iy,iz, p) = _BX(f, ix,b-iy-1,iz, p);
	  }
	  if (iz > -sw[2]) {
	    _BZ(f, ix,b+iy,iz, p) = _BZ(f, ix,b-iy-1,iz, p);
	  }
	}
	for (int iy = 1; iy < sw[1]; iy++) {
	  _BY(f, ix,my+iy,iz, p) = _BY(f, ix,my-iy-1,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_zh_open

static void
obndra_zh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       const int sw[3], const int ldims[3], int b, int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 2, p)) {
    int mz = ldims[2];
    for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iz = 0; iz < sw[2]; iz++) {
	  for (int m = 0; m < 5; m++) {
	    M3(f, m, ix,iy,b+iz, p) = M3(f, m, ix,iy,b-iz-1, p);
	  }
	  if (ix > -sw[0]) {
	    _BX(f, ix,iy,b+iz, p) = _BX(f, ix,iy,b-iz-1, p);
	  }
	  if (iy > -sw[1]) {
	    _BY(f, ix,iy,b+iz, p) = _BY(f, ix,iy,b-iz-1, p);
	  }
	}
	for (int iz = 1; iz < sw[2]; iz++) {
	  _BZ(f, ix,iy,mz+iz, p) = _BZ(f, ix,iy,mz-iz, p);
	}
      }
    }
  }
}

// ======================================================================
// cell-centered

// ----------------------------------------------------------------------	
// obndra_xl_open_cc
//
// set open boundary conditions at x-low for 5/10 moment fields

static void
obndra_xl_open_cc(struct ggcm_mhd_bnd *bnd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
    for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
      for (int ix = 0; ix > -sw[0]; ix--) {
	for (int m = 0; m < s_n_state; m++) {
	  M3(f, m, ix-1,iy,iz, p) = M3(f, m, ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_xh_open_cc

static void
obndra_xh_open_cc(struct ggcm_mhd *mhd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 0, p)) {
    int mx = ldims[0];
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
	for (int ix = 0; ix < sw[0]; ix++) {
	  for (int m = 0; m < s_n_state; m++) {
	    M3(f,m, mx+ix,iy,iz, p) = M3(f,m, mx-ix-1,iy,iz, p);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_yl_open_cc

static void
obndra_yl_open_cc(struct ggcm_mhd *mhd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 1, p)) {
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iy = 0; iy < sw[1]; iy++) {
	  for (int m = 0; m < s_n_state; m++) {
	    M3(f,m, ix,-1-iy,iz, p) = M3(f,m, ix,iy,iz, p);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_yh_open_cc

static void
obndra_yh_open_cc(struct ggcm_mhd *mhd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 1, p)) {
    int my = ldims[1];
    for (int iz = -sw[2]; iz < ldims[2] + sw[2]; iz++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iy = 0; iy < sw[1]; iy++) {
	  for (int m = 0; m < s_n_state; m++) {
	    M3(f,m, ix,my+iy,iz, p) = M3(f,m, ix,my-iy-1,iz, p);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_zl_open_cc

static void
obndra_zl_open_cc(struct ggcm_mhd *mhd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 2, p)) {
    for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iz = 0; iz < sw[2]; iz++) {
	  for (int m = 0; m < s_n_state; m++) {
	    M3(f,m, ix,iy,-1-iz, p) = M3(f,m, ix,iy,iz, p);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra_zh_open_cc

static void
obndra_zh_open_cc(struct ggcm_mhd *mhd, struct mrc_fld *f,
		  const int sw[3], const int ldims[3], int p)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 2, p)) {
    int mz = ldims[2];
    for (int iy = -sw[1]; iy < ldims[1] + sw[1]; iy++) {
      for (int ix = -sw[0]; ix < ldims[0] + sw[0]; ix++) {
	for (int iz = 0; iz < sw[2]; iz++) {
	  for (int m = 0; m < s_n_state; m++) {
	    M3(f,m, ix,iy,mz+iz, p) = M3(f,m, ix,iy,mz-iz-1, p);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndra
//
// set open fluid boundary conditions

static void
obndra(struct ggcm_mhd_bnd *bnd, struct mrc_fld *f, float bntim)
{
  struct ggcm_mhd_bnd_sub *sub = ggcm_mhd_bnd_sub(bnd);
  struct ggcm_mhd *mhd = bnd->mhd;

  const int *sw = mrc_fld_spatial_sw(f), *ldims = mrc_fld_spatial_dims(f);
  int bl[3], bh[3];
  for (int d = 0; d < 3; d++) {
    if (sub->do_legacy) {
      bl[d] = -1;
      bh[d] = ldims[d];
    } else {
      bl[d] = 0;
      bh[d] = ldims[d] - 1;
    }
  }
  
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    if (mrc_domain_at_boundary_lo(mhd->domain, 0, p)) {
      if (sub->apply_bndsw) {
	obndra_xl_bndsw(bnd, f, bntim, p);
      } else {
	if (MT_BGRID(MT) == MT_BGRID_CC) {
	  obndra_xl_open_cc(bnd, f, sw, ldims, p);
	} else {
	  assert(0);
	}
      }
    }

    // FIXME, the following two could be merged under
    // (MT_BGRID(MT) == MT_BGRID_CC)
    // if we're sure that the order doesn't matter
    if (MT == MT_GKEYLL) {
      obndra_yl_open_cc(mhd, f, sw, ldims, p);
      obndra_zl_open_cc(mhd, f, sw, ldims, p);
      
      obndra_xh_open_cc(mhd, f, sw, ldims, p);
      obndra_yh_open_cc(mhd, f, sw, ldims, p);
      obndra_zh_open_cc(mhd, f, sw, ldims, p);
    } else if (MT == MT_FCONS_CC) {
      obndra_xh_open_cc(mhd, f, sw, ldims, p);

      obndra_yl_open_cc(mhd, f, sw, ldims, p);
      obndra_yh_open_cc(mhd, f, sw, ldims, p);

      obndra_zl_open_cc(mhd, f, sw, ldims, p);
      obndra_zh_open_cc(mhd, f, sw, ldims, p);
    } else {
      obndra_xh_open(mhd, f, sw, ldims, bh[0], p);

      obndra_yl_open(mhd, f, sw, ldims, bl[1], p);
      obndra_yh_open(mhd, f, sw, ldims, bh[1], p);

      obndra_zl_open(mhd, f, sw, ldims, bl[2], p);
      obndra_zh_open(mhd, f, sw, ldims, bh[2], p);
    }
  }
}


// ----------------------------------------------------------------------
// obndrb_yl_open

static void
obndrb_yl_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       int m_t[3], int s_t[3], int l_n[3], int s_n[3], int p,
	       float *bdx3, float *bdy3, float *bdz3)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 1, p)) {
    for (int iz = -s_t[2]; iz < m_t[2] + s_t[2]; iz++) {
      for (int ix = -s_t[0]; ix < m_t[0] + s_t[0]; ix++) {
	for (int iy = l_n[1]; iy > l_n[1] - s_n[1]; iy--) {
	  _BY(f, ix,iy,iz, p) = BNDDIV_BY_L(ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndrb_zl_open

static void
obndrb_zl_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       int m_t[3], int s_t[3], int l_n[3], int s_n[3], int p,
	       float *bdx3, float *bdy3, float *bdz3)
{
  if (mrc_domain_at_boundary_lo(mhd->domain, 2, p)) {
    for (int iy = -s_t[1]; iy < m_t[1] + s_t[1]; iy++) {
      for (int ix = -s_t[0]; ix < m_t[0] + s_t[0]; ix++) {
	for (int iz = l_n[2]; iz > l_n[2] - s_n[2]; iz--) {
	  _BZ(f, ix,iy,iz, p) = BNDDIV_BZ_L(ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndrb_xh_open

static void
obndrb_xh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       int m_t[3], int s_t[3], int r_n[3], int s_n[3], int p,
	       float *bdx3, float *bdy3, float *bdz3)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 0, p)) {
    for (int iz = -s_t[2]; iz < m_t[2] + s_t[2]; iz++) {
      for (int iy = -s_t[1]; iy < m_t[1] + s_t[1]; iy++) {
	for (int ix = r_n[0]; ix < r_n[0] + s_n[0]; ix++) {
	  _BX(f, ix,iy,iz, p) = BNDDIV_BX_H(ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndrb_yh_open

static void
obndrb_yh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       int m_t[3], int s_t[3], int r_n[3], int s_n[3], int p,
	       float *bdx3, float *bdy3, float *bdz3)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 1, p)) {
    for (int iz = -s_t[2]; iz < m_t[2] + s_t[2]; iz++) {
      for (int ix = -s_t[0]; ix < m_t[0] + s_t[0]; ix++) {
	for (int iy = r_n[1]; iy < r_n[1] + s_n[1]; iy++) {
	  _BY(f, ix,iy,iz, p) = BNDDIV_BY_H(ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndrb_zh_open

static void
obndrb_zh_open(struct ggcm_mhd *mhd, struct mrc_fld *f,
	       int m_t[3], int s_t[3], int r_n[3], int s_n[3], int p,
	       float *bdx3, float *bdy3, float *bdz3)
{
  if (mrc_domain_at_boundary_hi(mhd->domain, 2, p)) {
    for (int iy = -s_t[1]; iy < m_t[1] + s_t[1]; iy++) {
      for (int ix = -s_t[0]; ix < m_t[0] + s_t[0]; ix++) {
	for (int iz = r_n[2]; iz < r_n[2] + s_n[2]; iz++) {
	  _BZ(f, ix,iy,iz, p) = BNDDIV_BZ_H(ix,iy,iz, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// obndrb

static void _mrc_unused
obndrb(struct ggcm_mhd_bnd *bnd, struct mrc_fld *f)
{
  struct ggcm_mhd_bnd_sub *sub = ggcm_mhd_bnd_sub(bnd);
  struct ggcm_mhd *mhd = bnd->mhd;

  const int *sw = mrc_fld_spatial_sw(f), *ldims = mrc_fld_spatial_dims(f);
#if MT == MT_FCONS_CC
  // tangential
  int m_t[3] = { ldims[0], ldims[1], ldims[2] }; // number of points
  int s_t[3] = { sw[0] - 1, sw[1] - 1, sw[2] - 1 }; // number of ghost points to fill
  // normal
  int l_n[3] = { -1, -1, -1 };
  int r_n[3] = { ldims[0], ldims[1], ldims[2] };
  int s_n[3] = { sw[0], sw[1], sw[2] };
#else
  // tangential
  int m_t[3] = { ldims[0], ldims[1], ldims[2] }; // number of points
  int s_t[3] = { sw[0] - 1, sw[1] - 1, sw[2] - 1 };
  // normal
  int l_n[3] = { 0, 0, 0 };
  int r_n[3] = { ldims[0], ldims[1], ldims[2] };
  int s_n[3] = { sw[0], sw[1], sw[2] }; // number of ghost points to fill
#endif

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    float *bdx3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
    float *bdy3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
    float *bdz3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);

    // assumes x1 bnd = fix, others = open
    if (!sub->do_legacy) {
      obndrb_yl_open(mhd, f, m_t, s_t, l_n, s_n, p, bdx3, bdy3, bdz3);
      obndrb_yh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
      
      obndrb_zl_open(mhd, f, m_t, s_t, l_n, s_n, p, bdx3, bdy3, bdz3);
      obndrb_zh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
      
      obndrb_xh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
    } else {
      obndrb_yl_open(mhd, f, m_t, s_t, l_n, s_n, p, bdx3, bdy3, bdz3);
      obndrb_zl_open(mhd, f, m_t, s_t, l_n, s_n, p, bdx3, bdy3, bdz3);

      obndrb_xh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
      obndrb_yh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
      obndrb_zh_open(mhd, f, m_t, s_t, r_n, s_n, p, bdx3, bdy3, bdz3);
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sub_fill_ghosts

static void
ggcm_mhd_bnd_sub_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
			     float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;
  assert(mhd);

  static int PR;
  if (!PR) {
    PR = prof_register(__FUNCTION__, 1., 0, 0);
  }
  prof_start(PR);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);
  assert(mhd_type == MT);
 
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  obndra(bnd, f, bntim);
#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
  assert(f->_aos && f->_c_order);
#else
  obndrb(bnd, f);
#endif
  mrc_fld_put_as(f, fld);

  prof_stop(PR);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sub_setup

static void
ggcm_mhd_bnd_sub_setup(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));

  if (!mhd->bnd_mask) {
    return;
  }

  struct mrc_fld *bnd_mask = mrc_fld_get_as(mhd->bnd_mask, FLD_TYPE);
  const int *dims = mrc_fld_spatial_dims(bnd_mask);

  for (int p = 0; p < mrc_fld_nr_patches(bnd_mask); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    mrc_fld_foreach(bnd_mask, ix,iy,iz, 2, 2) {
      int i[3] = { ix, iy, iz };
      for (int d = 0; d < 3; d++) {
	if (mrc_domain_at_boundary_lo(mhd->domain, d, p)) {
	  if (i[d] < 0) {
	    M3(bnd_mask, 0, ix,iy,iz, p) = 1.f;
	  } else if (i[d] == 0) {
	    M3(bnd_mask, 0, ix,iy,iz, p) = 2.f;
	  }
	}
 
	if (mrc_domain_at_boundary_hi(mhd->domain, d, p)) {
	  if (i[d] >= dims[d]) {
	    M3(bnd_mask, 0, ix,iy,iz, p) = 1.f;
	  } else if (i[d] == dims[d] - 1) {
	    M3(bnd_mask, 0, ix,iy,iz, p) = 2.f;
	  }
	}
      }
   } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(bnd_mask, mhd->bnd_mask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sub_destroy

static void
ggcm_mhd_bnd_sub_destroy(struct ggcm_mhd_bnd *bnd)
{
  pde_free();
}


// ----------------------------------------------------------------------
// ggcm_mhd_bnd inflow description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bnd_sub, x)
static struct param ggcm_mhd_bnd_sub_descr[] = {
  { "rr"           , VAR(bnvals[RR])        , PARAM_DOUBLE(1.) },
  { "pp"           , VAR(bnvals[PP])        , PARAM_DOUBLE(1.) },
  { "vx"           , VAR(bnvals[VX])        , PARAM_DOUBLE(0.) },
  { "vy"           , VAR(bnvals[VY])        , PARAM_DOUBLE(0.) },
  { "vz"           , VAR(bnvals[VZ])        , PARAM_DOUBLE(0.) },
  { "bx"           , VAR(bnvals[BX])        , PARAM_DOUBLE(0.) },
  { "by"           , VAR(bnvals[BY])        , PARAM_DOUBLE(0.) },
  { "bz"           , VAR(bnvals[BZ])        , PARAM_DOUBLE(0.) },
  { "apply_bndsw"  , VAR(apply_bndsw)       , PARAM_BOOL(true) },
  { "do_legacy"    , VAR(do_legacy)         , PARAM_BOOL(false)},

  {},
};
#undef VAR

// ======================================================================
// ggcm_mhd_bnd subclass "inoutflow"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_inoutflow = {
  .name             = ggcm_mhd_bnd_sub_name,
  .size             = sizeof(struct ggcm_mhd_bnd_sub),
  .param_descr      = ggcm_mhd_bnd_sub_descr,
  .setup            = ggcm_mhd_bnd_sub_setup,
  .destroy          = ggcm_mhd_bnd_sub_destroy,
  .fill_ghosts      = ggcm_mhd_bnd_sub_fill_ghosts,
};

