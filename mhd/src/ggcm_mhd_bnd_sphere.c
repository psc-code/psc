
#include "ggcm_mhd_bnd_private.h"

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>

#include <mrc_fld_as_float.h>
#include <mrc_domain.h>
#include <mrc_bits.h>
#include <math.h>

#define MT MT_FULLY_CONSERVATIVE

enum {
  FIXED_RR,
  FIXED_PP,
  FIXED_VX,
  FIXED_VY,
  FIXED_VZ,
  FIXED_BX,
  FIXED_BY,
  FIXED_BZ,
  FIXED_NR,
};


// FIXME, consolidate with ggcm_mhd_iono

// ======================================================================
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_sphere {
  // params
  double radius;

  // state
  double min_dr;
  double r1;
  double r2;

  // maps
  // for managing cell-centered ghost points
  int cc_n_map;
  struct mrc_fld *cc_mhd_imap;  // ghost cell # -> (ix,iy,iz,p)

  // constant values to set
  double bnvals[FIXED_NR];
};

#define ggcm_mhd_bnd_sphere(bnd) mrc_to_subobj(bnd, struct ggcm_mhd_bnd_sphere)

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_map_find_dr
//
// find minimum cell size (over all of the domain -- FIXME?)

static void
ggcm_mhd_bnd_map_find_dr(struct ggcm_mhd_bnd *bnd, double *dr)
{
  // FIXME, it'd make sense to base this on the ionosphere boundary region
  // (e.g., box of +/- 7 RE in all directions, as later used in determining
  // r1, r2). It shouldn't really hurt if the dr determined here is too small,
  // though it'll slow down finding the proper r1, r2.

  struct ggcm_mhd *mhd = bnd->mhd;

  double min_dr = 1.e30;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    for (int d = 0; d < 3; d++) {
      for (int i = 0; i < info.ldims[d]; i++) {
	min_dr = fmin(min_dr, MRC_MCRD(crds, d, i, p) - MRC_MCRD(crds, d, i-1, p));
      }
    }
  }
  MPI_Allreduce(&min_dr, dr, 1, MPI_DOUBLE, MPI_MIN, ggcm_mhd_bnd_comm(bnd));
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_map_find_r1_r2
//
// determines r1 to be small enough so that
// +/- 2 grid points around each cell with center outside of r2
// are outside (ie., their centers) the smaller r1 sphere

static void
ggcm_mhd_bnd_map_find_r1_r2(struct ggcm_mhd_bnd *bnd,
			    double *p_r1, double *p_r2)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);

  struct ggcm_mhd *mhd = bnd->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double dr = sub->min_dr;
  double r2 = sub->radius;
  double r1 = r2 - dr;

 loop2:
  // check if r1, r2 are ok
  for (int p = 0; p < mrc_domain_nr_patches(mhd->domain); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    for (int iz = 0; iz < info.ldims[2]; iz++) {
      for (int iy = 0; iy < info.ldims[1]; iy++) {
	for (int ix = 0; ix < info.ldims[0]; ix++) {
	  double xx = MRC_MCRDX(crds, ix, p);
	  double yy = MRC_MCRDY(crds, iy, p);
	  double zz = MRC_MCRDZ(crds, iz, p);
	  
	  double rr = sqrt(sqr(xx) + sqr(yy) + sqr(zz));
	  if (rr <= r2) continue;
	  
	  for (int jz = -2; jz <= 2; jz++) {
	    for (int jy = -2; jy <= 2; jy++) {
	      for (int jx = -2; jx <= 2; jx++) {
		double xxx = MRC_MCRDX(crds, ix+jx, p);
		double yyy = MRC_MCRDY(crds, iy+jy, p);
		double zzz = MRC_MCRDZ(crds, iz+jz, p);
		double rrr = sqrt(sqr(xxx) + sqr(yyy) + sqr(zzz));
		if (rrr < r1) {
		  r1 -= .01; // FIXME, hardcoded number
		  goto loop2;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  MPI_Allreduce(&r1, p_r1, 1, MPI_DOUBLE, MPI_MIN, ggcm_mhd_bnd_comm(bnd));

  *p_r2 = r2;
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_map_find_cc_n_map

static void
ggcm_mhd_bnd_map_find_cc_n_map(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd *mhd = bnd->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double r1 = sub->r1, r2 = sub->r2;
  assert(r1 > 0.);

  int cc_n_map = 0;
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    int gdims[3];
    mrc_domain_get_global_dims(mhd->domain, gdims);
    // cell-centered
    int sw[3] = { 2, 2, 2 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	sw[d] = 0;
      }
    }
    for (int jz = -sw[2]; jz < info.ldims[2] + sw[2]; jz++) {
      for (int jy = -sw[1]; jy < info.ldims[1] + sw[1]; jy++) {
	for (int jx = -sw[0]; jx < info.ldims[0] + sw[0]; jx++) {
	  float xx = MRC_MCRDX(crds, jx, p);
	  float yy = MRC_MCRDY(crds, jy, p);
	  float zz = MRC_MCRDZ(crds, jz, p);
	  float rr = sqrtf(sqr(xx) + sqr(yy) + sqr(zz));
	  if (rr < r1 || rr > r2) continue;
	  cc_n_map++;
	}
      }
    }
  }
  sub->cc_n_map = cc_n_map;
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_map_cc

static void
ggcm_mhd_bnd_map_cc(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);

  struct ggcm_mhd *mhd = bnd->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double r1 = sub->r1, r2 = sub->r2;

  // compute e-field mapping coefficients

  int cc_n_map = 0;
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    int gdims[3];
    mrc_domain_get_global_dims(mhd->domain, gdims);
    // cell-centered
    int sw[3] = { 2, 2, 2 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	sw[d] = 0;
      }
    }
    for (int jz = -sw[2]; jz < info.ldims[2] + sw[2]; jz++) {
      for (int jy = -sw[1]; jy < info.ldims[1] + sw[1]; jy++) {
	for (int jx = -sw[0]; jx < info.ldims[0] + sw[0]; jx++) {
	  double xx = MRC_MCRDX(crds, jx, p);
	  double yy = MRC_MCRDY(crds, jy, p);
	  double zz = MRC_MCRDZ(crds, jz, p);
	  double rr = sqrtf(sqr(xx) + sqr(yy) + sqr(zz));
	  if (rr < r1 || rr > r2) continue;
	  
	  MRC_I2(sub->cc_mhd_imap, 0, cc_n_map) = jx;
	  MRC_I2(sub->cc_mhd_imap, 1, cc_n_map) = jy;
	  MRC_I2(sub->cc_mhd_imap, 2, cc_n_map) = jz;
	  MRC_I2(sub->cc_mhd_imap, 3, cc_n_map) = p;

	  cc_n_map++;
	}
      }
    }
  }

  assert(cc_n_map == sub->cc_n_map);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_setup_flds

static void
ggcm_mhd_bnd_sphere_setup_flds(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);

  ggcm_mhd_bnd_map_find_cc_n_map(bnd);
  mprintf("cc_n_map %d\n", sub->cc_n_map);

  // cell-centered

  mrc_fld_set_type(sub->cc_mhd_imap, "int");
  mrc_fld_set_param_int_array(sub->cc_mhd_imap, "dims", 2, (int[2]) { 4, sub->cc_n_map });
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_setup

static void
ggcm_mhd_bnd_sphere_setup(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);

  ggcm_mhd_bnd_map_find_dr(bnd, &sub->min_dr);
  mprintf("min_dr %g\n", sub->min_dr);
  ggcm_mhd_bnd_map_find_r1_r2(bnd, &sub->r1, &sub->r2);
  mprintf("r1 %g r2 %g\n", sub->r1, sub->r2);
  ggcm_mhd_bnd_sphere_setup_flds(bnd);
  ggcm_mhd_bnd_setup_member_objs_sub(bnd);
  ggcm_mhd_bnd_map_cc(bnd);
}

// ----------------------------------------------------------------------
// sphere_fill_ghosts_mhd_do

static void
sphere_fill_ghosts_mhd_do(struct mrc_fld *fld,
    int cc_n_map, struct mrc_fld*cc_mhd_imap,
    double bnvals[FIXED_NR], int m, float bntim, float gamm)
{
  double rvx = bnvals[FIXED_RR] * bnvals[FIXED_VX];
  double rvy = bnvals[FIXED_RR] * bnvals[FIXED_VY];
  double rvz = bnvals[FIXED_RR] * bnvals[FIXED_VZ];

  double vvbn  = sqr(bnvals[FIXED_VX]) + sqr(bnvals[FIXED_VY]) + sqr(bnvals[FIXED_VZ]);
  double uubn  = .5f * (bnvals[FIXED_RR]*vvbn) + bnvals[FIXED_PP] / (gamm - 1.f);
  double b2bn  = sqr(bnvals[FIXED_BX]) + sqr(bnvals[FIXED_BY]) + sqr(bnvals[FIXED_BZ]);
  double eebn = uubn + .5 * b2bn;

  for (int i = 0; i < cc_n_map; i++) {
    int ix = MRC_I2(cc_mhd_imap, 0, i);
    int iy = MRC_I2(cc_mhd_imap, 1, i);
    int iz = MRC_I2(cc_mhd_imap, 2, i);
    int p  = MRC_I2(cc_mhd_imap, 3, i);

    M3 (fld, m + RR,  ix,iy,iz, p) = bnvals[FIXED_RR];
    M3 (fld, m + RVX, ix,iy,iz, p) = rvx;
    M3 (fld, m + RVY, ix,iy,iz, p) = rvy;
    M3 (fld, m + RVZ, ix,iy,iz, p) = rvz;
    if (MT == MT_SEMI_CONSERVATIVE ||
        MT == MT_SEMI_CONSERVATIVE_GGCM) {
      M3(fld, m + UU , ix,iy,iz, p) = uubn;
    } else if (MT == MT_FULLY_CONSERVATIVE) {
      M3(fld, m + EE , ix,iy,iz, p) = eebn;
    } else {
      assert(0);
    }
    M3(fld, m + BX , ix,iy,iz, p) = bnvals[FIXED_BX];
    M3(fld, m + BY , ix,iy,iz, p) = bnvals[FIXED_BY];
    M3(fld, m + BZ , ix,iy,iz, p) = bnvals[FIXED_BZ];
  }
}


// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_fill_ghosts

static void
ggcm_mhd_bnd_sphere_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld_base,
			      int m, float bntim)
{
  struct ggcm_mhd_bnd_sphere *sub = ggcm_mhd_bnd_sphere(bnd);
  struct ggcm_mhd *mhd = bnd->mhd;

  if (sub->cc_n_map == 0) {
    return;
  }

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  sphere_fill_ghosts_mhd_do(fld, sub->cc_n_map, sub->cc_mhd_imap,
      sub->bnvals, m, bntim, mhd->par.gamm);

  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd "sphere" subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_bnd_sphere, x)
static struct param ggcm_mhd_bnd_sphere_descr[] = {
  { "radius"          , VAR(radius)          , PARAM_DOUBLE(1.)          },

  { "min_dr"          , VAR(min_dr)          , MRC_VAR_FLOAT             },
  { "r1"              , VAR(r1)              , MRC_VAR_FLOAT             },
  { "r2"              , VAR(r2)              , MRC_VAR_FLOAT             },
  { "cc_n_map"        , VAR(cc_n_map)        , MRC_VAR_INT               },

  { "cc_mhd_imap"     , VAR(cc_mhd_imap)     , MRC_VAR_OBJ(mrc_fld)      },

  { "rr"              , VAR(bnvals[FIXED_RR]), PARAM_DOUBLE(1.) },
  { "pp"              , VAR(bnvals[FIXED_PP]), PARAM_DOUBLE(1.) },
  { "vx"              , VAR(bnvals[FIXED_VX]), PARAM_DOUBLE(0.) },
  { "vy"              , VAR(bnvals[FIXED_VY]), PARAM_DOUBLE(0.) },
  { "vz"              , VAR(bnvals[FIXED_VZ]), PARAM_DOUBLE(0.) },
  { "bx"              , VAR(bnvals[FIXED_BX]), PARAM_DOUBLE(0.) },
  { "by"              , VAR(bnvals[FIXED_BY]), PARAM_DOUBLE(0.) },
  { "bz"              , VAR(bnvals[FIXED_BZ]), PARAM_DOUBLE(0.) },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_bnd subclass "sphere"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_sphere = {
  .name             = "sphere",
  .size             = sizeof(struct ggcm_mhd_bnd_sphere),
  .param_descr      = ggcm_mhd_bnd_sphere_descr,
  .setup            = ggcm_mhd_bnd_sphere_setup,
  .fill_ghosts      = ggcm_mhd_bnd_sphere_fill_ghosts,
};

