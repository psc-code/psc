
#include "ggcm_mhd_bndsw.h"
#include "ggcm_mhd_crds.h"
#include "mrc_domain.h"

// FIXME
// The B boundary conditions have numerous issues:
//
// - Originally (and still in the Fortran version), there are issues
//   in actually getting div B = 0 as desired.  
//
// - Originally (and still in the Fortran version), there is an asymmetry in
//   the lower vs upper boundaries, in particular in y, z.
//   The code would set, e.g., By[-2], By[-1] and By[my], By[my+1].
//   In the original staggering, ie., the By on the boundary and one cell to
//   the left would be set on the lower side, but on the upper side, the code would
//   set the two By values to the right of the boundary, not the one on the wall
//
// - With the new staggering, we could only do the opposite, which is still
//   asymmetric, and also would imply different b.c. depending on the
//   the staggering chosen. So instead, the ggcm_mhd_bnd "c" and "c2" now
//   set the values right on the boundary and the next one outside of the
//   boundary, which is different from Fortran, but symmetric and
//   allows us to get the same results independently of which staggering is
//   chosen (well, almost, it seems there's another rather small error
//   appearing after a couple of steps, which may or may not be real)
//
// - It's not clear how many ghost points we need in the first place, and how
//   this all interacts with the bpush update (I think bpush updates a boundary
//   point, which however afterwards will be fixed up here, so should be okay)
//   The current implementation of eta J in Ohm's Law probably really needs two
//   ghost points.
//
// - The "Set div B = 0" implementation at corners is iffy (order dependent), it
//   will, e.g., set By based on Bx at the wall (which is not well defined), and
//   then Bx at the wall based on that By. It's unlikely to really matter, though.
//   (But it might be the reason for the remaining small discrepancies between the
//   two staggerings.)

#if SHIFT == -1

// DIFF bdy2[] index etc was inconsistent, now actually divg free
// FIXME, which of these ghost points are actually used? / loop limits

#define BNDDIV_BY_L(ix, iy, iz, p)					\
  (M3(f, mm+1, ix,iy+1,iz, p) + (1. / bdy3[iy+1]) *			\
   (bdx3[ix] * (M3(f, mm+0, ix  ,iy+1,iz  , p) -			\
		M3(f, mm+0, ix-1,iy+1,iz  , p)) +			\
    bdz3[iz] * (M3(f, mm+2, ix  ,iy+1,iz  , p) -			\
		M3(f, mm+2, ix  ,iy+1,iz-1, p))))
#define BNDDIV_BZ_L(ix, iy, iz, p)					\
  (M3(f, mm+2, ix,iy,iz+1, p) + (1. / bdz3[iz+1]) *			\
   (bdx3[ix] * (M3(f, mm+0, ix  ,iy  ,iz+1, p) -			\
		M3(f, mm+0, ix-1,iy  ,iz+1, p)) +			\
    bdy3[iy] * (M3(f, mm+1, ix  ,iy  ,iz+1, p) -			\
		M3(f, mm+1, ix  ,iy-1,iz+1, p))))

#define BNDDIV_BX_H(ix, iy, iz, p)					\
  (M3(f, mm+0, ix-1,iy,iz, p) - (1. / bdx3[ix]) *			\
   (bdy3[iy] * (M3(f, mm+1, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+1, ix  ,iy-1,iz  , p)) +			\
    bdz3[iz] * (M3(f, mm+2, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+2, ix  ,iy  ,iz-1, p))))
#define BNDDIV_BY_H(ix, iy, iz, p)					\
  (M3(f, mm+1, ix,iy-1,iz, p) - (1. / bdy3[iy]) *			\
   (bdx3[ix] * (M3(f, mm+0, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+0, ix-1,iy  ,iz  , p)) +			\
    bdz3[iz] * (M3(f, mm+2, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+2, ix  ,iy  ,iz-1, p))))
#define BNDDIV_BZ_H(ix, iy, iz, p)					\
  (M3(f, mm+2, ix,iy,iz-1, p) - (1. / bdz3[iz]) *			\
   (bdx3[ix] * (M3(f, mm+0, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+0, ix-1,iy  ,iz  , p)) +			\
    bdy3[iy] * (M3(f, mm+1, ix  ,iy  ,iz  , p) -			\
		M3(f, mm+1, ix  ,iy-1,iz  , p))))

#elif SHIFT == 0

// DIFF bdy2[] index etc was inconsistent, now actually divg free
// FIXME, which of these ghost points are actually used? / loop limits

#define BNDDIV_BY_L(ix, iy, iz, p)					\
  (M3(f, mm+1, ix,iy+1,iz, p) +						\
   bdx3[ix]/bdy3[iy] * (M3(f, mm+0, ix+1,iy,iz  , p) -			\
			M3(f, mm+0, ix  ,iy,iz  , p)) +			\
   bdz3[iz]/bdy3[iy] * (M3(f, mm+2, ix  ,iy,iz+1, p) -			\
			M3(f, mm+2, ix  ,iy,iz  , p)))
#define BNDDIV_BZ_L(ix, iy, iz, p)					\
  (M3(f, mm+2, ix,iy,iz+1, p) +						\
   bdx3[ix]/bdz3[iz] * (M3(f, mm+0, ix+1,iy  ,iz, p) -			\
			M3(f, mm+0, ix  ,iy  ,iz, p)) +			\
   bdy3[iy]/bdz3[iz] * (M3(f, mm+1, ix  ,iy+1,iz, p) -			\
			M3(f, mm+1, ix  ,iy  ,iz, p)))

#define BNDDIV_BX_H(ix, iy, iz, p)					\
  (M3(f, mm+0, ix-1,iy,iz, p) -						\
   bdy3[iy]/bdx3[ix-1] * (M3(f, mm+1, ix-1,iy+1,iz  , p) -		\
			  M3(f, mm+1, ix-1,iy  ,iz  , p)) -		\
   bdz3[iz]/bdx3[ix-1] * (M3(f, mm+2, ix-1,iy  ,iz+1, p) -		\
			  M3(f, mm+2, ix-1,iy  ,iz  , p)))
#define BNDDIV_BY_H(ix, iy, iz, p)					\
  (M3(f, mm+1, ix,iy-1,iz, p) -						\
   bdx3[ix]/bdy3[iy-1] * (M3(f, mm+0, ix+1,iy-1,iz  , p) -		\
			  M3(f, mm+0, ix  ,iy-1,iz  , p)) -		\
   bdz3[iz]/bdy3[iy-1] * (M3(f, mm+2, ix  ,iy-1,iz+1, p) -		\
			  M3(f, mm+2, ix  ,iy-1,iz  , p)))
#define BNDDIV_BZ_H(ix, iy, iz, p)					\
  (M3(f, mm+2, ix,iy,iz-1, p) -						\
   bdx3[ix]/bdz3[iz-1] * (M3(f, mm+0, ix+1,iy  ,iz-1, p) -		\
			  M3(f, mm+0, ix  ,iy  ,iz-1, p)) -		\
   bdy3[iy]/bdz3[iz-1] * (M3(f, mm+1, ix  ,iy+1,iz-1, p) -		\
			  M3(f, mm+1, ix  ,iy  ,iz-1, p)))

#else

#error unknown SHIFT

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
bnd_sw(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p, double vals[8], float bntim)
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  float xx[3] = { MRC_MCRDX(crds, ix, p),
		  MRC_MCRDY(crds, iy, p),
		  MRC_MCRDZ(crds, iz, p), };
  float bn[SW_NR];

  static struct ggcm_mhd_bndsw *bndsw;
  if (!bndsw) {
    bndsw = ggcm_mhd_get_var_obj(mhd, "bndsw");
  }
  ggcm_mhd_bndsw_at(bndsw, bntim, xx, bn);

  float vvbn  = sqr(bn[SW_VX]) + sqr(bn[SW_VY]) + sqr(bn[SW_VZ]);
  float uubn  = .5f * (bn[SW_RR]*vvbn) + bn[SW_PP] / (mhd->par.gamm - 1.f);

  vals[RR ] = bn[SW_RR];
  vals[RVX] = bn[SW_RR] * bn[SW_VX];
  vals[RVY] = bn[SW_RR] * bn[SW_VY];
  vals[RVZ] = bn[SW_RR] * bn[SW_VZ];
  vals[UU ] = uubn;
  vals[BX ] = bn[SW_BX];
  vals[BY ] = bn[SW_BY];
  vals[BZ ] = bn[SW_BZ];
}

// ----------------------------------------------------------------------	
// obndra
//
// set fluid boundary conditions at inflow boundary

static void
obndra(struct ggcm_mhd *mhd, struct mrc_fld *f, int mm, float bntim)
{
  static int PR;
  if (!PR) {
    PR = prof_register(__FUNCTION__, 1., 0, 0);
  }
  prof_start(PR);
  int sw = f->_nr_ghosts;
  int mx = mhd->im[0], my = mhd->im[1], mz = mhd->im[2];

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    if (mrc_domain_at_boundary_lo(mhd->domain, 0, p)) {
      for (int iz = -sw; iz < mz + sw; iz++) {
	for (int iy = -sw; iy < my + sw; iy++) {
	  for (int ix = -sw; ix < 0; ix++) {
	    double vals[8];
	    bnd_sw(mhd, ix, iy, iz, p, vals, bntim);
	    if (MT == MT_FULLY_CONSERVATIVE) {
	      mrc_fld_data_t b2  = sqr(vals[BX]) + sqr(vals[BY]) + sqr(vals[BZ]);
 	      vals[EE] += .5 * b2;
	    }
	    for (int m = 0; m < 8; m++) {
	      M3(f, mm + m, ix,iy,iz, p) = vals[m];
	    }
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_lo(mhd->domain, 1, p)) {
      for (int iz = -sw; iz < mz + sw; iz++) {
	for (int ix = -sw; ix < mx + sw; ix++) {
	  for (int iy = 0; iy > -sw; iy--) {
	    for (int m = mm; m < mm + 5; m++) {
	      M3(f,m, ix,iy-1,iz, p) = M3(f,m, ix,iy,iz, p);
	    }
	    M3(f,mm + BX, ix,iy-1,iz, p) = M3(f,mm + BX, ix,iy,iz, p);
	    M3(f,mm + BZ, ix,iy-1,iz, p) = M3(f,mm + BZ, ix,iy,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_lo(mhd->domain, 2, p)) {
      for (int iy = -sw; iy < my + sw; iy++) {
	for (int ix = -sw; ix < mx + sw; ix++) {
	  for (int iz = 0; iz > -sw; iz--) {
	    for (int m = mm; m < mm + 5; m++) {
	      M3(f,m, ix,iy,iz-1, p) = M3(f,m, ix,iy,iz, p);
	    }
	    M3(f,mm + BX, ix,iy,iz-1, p) = M3(f,mm + BX, ix,iy,iz, p);
	    M3(f,mm + BY, ix,iy,iz-1, p) = M3(f,mm + BY, ix,iy,iz, p);
	  }
	}
      }
    }
    
    if (mrc_domain_at_boundary_hi(mhd->domain, 0, p)) {
      for (int iz = -sw; iz < mz + sw; iz++) {
	for (int iy = -sw; iy < my + sw; iy++) {
	  for (int ix = mx; ix < mx + sw; ix++) {
	    for (int m = mm; m < mm + 5; m++) {
	      M3(f,m, ix,iy,iz, p) = M3(f,m, ix-1,iy,iz, p);
	    }
	    M3(f,mm + BY, ix,iy,iz, p) = M3(f,mm + BY, ix-1,iy,iz, p);
	    M3(f,mm + BZ, ix,iy,iz, p) = M3(f,mm + BZ, ix-1,iy,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_hi(mhd->domain, 1, p)) {
      for (int iz = -sw; iz < mz + sw; iz++) {
	for (int ix = -sw; ix < mx + sw; ix++) {
	  for (int iy = my; iy < my + sw; iy++) {
	    for (int m = mm; m < mm + 5; m++) {
	      M3(f,m, ix,iy,iz, p) = M3(f,m, ix,iy-1,iz, p);
	    }
	    M3(f,mm + BX, ix,iy,iz, p) = M3(f,mm + BX, ix,iy-1,iz, p);
	    M3(f,mm + BZ, ix,iy,iz, p) = M3(f,mm + BZ, ix,iy-1,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_hi(mhd->domain, 2, p)) {
      for (int iy = -sw; iy < my + sw; iy++) {
	for (int ix = -sw; ix < mx + sw; ix++) {
	  for (int iz = mz; iz < mz + sw; iz++) {
	    for (int m = mm; m < mm + 5; m++) {
	      M3(f,m, ix,iy,iz, p) = M3(f,m, ix,iy,iz-1, p);
	    }
	    M3(f,mm + BX, ix,iy,iz, p) = M3(f,mm + BX, ix,iy,iz-1, p);
	    M3(f,mm + BY, ix,iy,iz, p) = M3(f,mm + BY, ix,iy,iz-1, p);
	  }
	}
      }
    }
  }
  prof_stop(PR);
}

static void
obndrb(struct ggcm_mhd *mhd, struct mrc_fld *f, int mm)
{
  static int PR;
  if (!PR) {
    PR = prof_register(__FUNCTION__, 1., 0, 0);
  }
  prof_start(PR);

  int sw = f->_nr_ghosts;
  int mx = mhd->im[0], my = mhd->im[1], mz = mhd->im[2];

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    float *bdx3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
    float *bdy3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
    float *bdz3 = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);

    // assumes x1 bnd = fix, others = open
    if (mrc_domain_at_boundary_lo(mhd->domain, 1, p)) {
      for (int iz = -sw - SHIFT; iz < mz + sw - 1 - SHIFT; iz++) {
	for (int ix = -sw - SHIFT; ix < mx + sw - 1 - SHIFT; ix++) {
	  for (int iy = 0; iy > -sw; iy--) {
	    M3(f, mm+1, ix,SHIFT + iy,iz, p) = BNDDIV_BY_L(ix,SHIFT + iy,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_lo(mhd->domain, 2, p)) {
      for (int iy = -sw - SHIFT; iy < my + sw - 1 - SHIFT; iy++) {
	for (int ix = -sw - SHIFT; ix < mx + sw - 1 - SHIFT; ix++) {
	  for (int iz = 0; iz > -sw; iz--) {
	    M3(f, mm+2, ix,iy,SHIFT + iz, p) = BNDDIV_BZ_L(ix,iy,SHIFT + iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_hi(mhd->domain, 0, p)) {
      for (int iz = -sw - SHIFT; iz < mz + sw - 1 - SHIFT; iz++) {
	for (int iy = -sw - SHIFT; iy < my + sw - 1 - SHIFT; iy++) {
	  for (int ix = mx; ix < mx + sw; ix++) {
	    M3(f, mm+0, SHIFT + ix,iy,iz, p) = BNDDIV_BX_H(SHIFT + ix,iy,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_hi(mhd->domain, 1, p)) {
      for (int iz = -sw - SHIFT; iz < mz + sw - 1 - SHIFT; iz++) {
	for (int ix = -sw - SHIFT; ix < mx + sw - 1 - SHIFT; ix++) {
	  for (int iy = my; iy < my + sw; iy++) {
	    M3(f, mm+1, ix,SHIFT + iy,iz, p) = BNDDIV_BY_H(ix,SHIFT + iy,iz, p);
	  }
	}
      }
    }
    if (mrc_domain_at_boundary_hi(mhd->domain, 2, p)) {
      for (int iy = -sw - SHIFT; iy < my + sw - 1 - SHIFT; iy++) {
	for (int ix = -sw - SHIFT; ix < mx + sw - 1 - SHIFT; ix++) {
	  for (int iz = mz; iz < mz + sw; iz++) {
	    M3(f, mm+2, ix,iy,SHIFT + iz, p) = BNDDIV_BZ_H(ix,iy,SHIFT + iz, p);
	  }
	}
      }
    }
  }
  prof_stop(PR);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sub_fill_ghosts

static void
ggcm_mhd_bnd_sub_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
			     int m, float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;
  assert(mhd);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);
  assert(mhd_type == MT);

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  assert(m == 0 || m == 8);
  obndra(mhd, f, m, bntim);
  obndrb(mhd, f, m + BX);
  mrc_fld_put_as(f, fld);
}

