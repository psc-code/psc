
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_bits.h>
#include <mrc_fld_as_double.h>

#include <math.h>

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_is_bnd

static bool
ggcm_mhd_bnd_sphere_map_is_bnd(struct ggcm_mhd_bnd_sphere_map *map, double xx[3])
{
  return (sqr(xx[0]) + sqr(xx[1]) + sqr(xx[2])) < sqr(map->radius);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_is_ghost

static bool
ggcm_mhd_bnd_sphere_map_is_ghost(struct ggcm_mhd_bnd_sphere_map *map, double xx[3])
{
  double rr = sqr(xx[0]) + sqr(xx[1]) + sqr(xx[2]);
  return (rr < sqr(map->radius) && rr > sqr(map->r1));
    
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_set_ymask

static void
ggcm_mhd_bnd_sphere_map_set_ymask(struct ggcm_mhd_bnd_sphere_map *map)
{
  assert(map->mhd->ymask);
  struct mrc_fld *ymask = mrc_fld_get_as(map->mhd->ymask, FLD_TYPE);

  struct mrc_crds *crds = mrc_domain_get_crds(map->mhd->domain);  

  for (int p = 0; p < mrc_fld_nr_patches(ymask); p++) {
    mrc_fld_foreach(ymask, ix,iy,iz, 2, 2) {
      double xx[3] = { MRC_MCRDX(crds, ix, p),
		       MRC_MCRDY(crds, iy, p),
		       MRC_MCRDZ(crds, iz, p) };
      if (ggcm_mhd_bnd_sphere_map_is_bnd(map, xx)) {
	M3(ymask, 0, ix,iy,iz, p) = 0.;
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(ymask, map->mhd->ymask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_set_bnd_mask

static void
ggcm_mhd_bnd_sphere_map_set_bnd_mask(struct ggcm_mhd_bnd_sphere_map *map)
{
  assert(map->mhd->bnd_mask);
  struct mrc_fld *bnd_mask = mrc_fld_get_as(map->mhd->bnd_mask, FLD_TYPE);

  struct mrc_crds *crds = mrc_domain_get_crds(map->mhd->domain);

  for (int p = 0; p < mrc_fld_nr_patches(bnd_mask); p++) {
    mrc_fld_foreach(bnd_mask, ix,iy,iz, 2, 2) {
      double xx[3] = { MRC_MCRDX(crds, ix, p),
		       MRC_MCRDY(crds, iy, p),
		       MRC_MCRDZ(crds, iz, p) };
      if (ggcm_mhd_bnd_sphere_map_is_bnd(map, xx)) {
	M3(bnd_mask, 0, ix,iy,iz, p) = 1.;
      }
    } mrc_fld_foreach_end;
  }
  
  mrc_fld_put_as(bnd_mask, map->mhd->bnd_mask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_dr
//
// find minimum cell size (over all of the domain -- FIXME?)

static void
ggcm_mhd_bnd_sphere_map_find_dr(struct ggcm_mhd_bnd_sphere_map *map, double *dr)
{
  // FIXME, it'd make sense to base this on the ionosphere boundary region
  // (e.g., box of +/- 7 RE in all directions, as later used in determining
  // r1, r2). It shouldn't really hurt if the dr determined here is too small,
  // though it'll slow down finding the proper r1, r2.

  struct ggcm_mhd *mhd = map->mhd;

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
  MPI_Allreduce(&min_dr, dr, 1, MPI_DOUBLE, MPI_MIN, ggcm_mhd_comm(mhd));
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_r1
//
// determines r1 to be small enough so that
// +/- 2 grid points around each cell with center outside of r2
// are outside (ie., their centers) the smaller r1 sphere

static void
ggcm_mhd_bnd_sphere_map_find_r1(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double dr = map->min_dr;
  double r2 = map->radius;
  double r1 = r2 - dr;

 loop2:
  // check if r1 is ok
  for (int p = 0; p < mrc_domain_nr_patches(mhd->domain); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    for (int iz = 0; iz < info.ldims[2]; iz++) {
      for (int iy = 0; iy < info.ldims[1]; iy++) {
	for (int ix = 0; ix < info.ldims[0]; ix++) {
	  double xx[3] = { MRC_MCRDX(crds, ix, p),
			   MRC_MCRDY(crds, iy, p),
			   MRC_MCRDZ(crds, iz, p), };
	  if (ggcm_mhd_bnd_sphere_map_is_bnd(map, xx)) {
	    // if we're inside the sphere, nothing to check
	    continue;
	  }
	  double rr = sqrt(sqr(xx[0]) + sqr(xx[1]) + sqr(xx[2]));
	  if (rr > 2. * r2) {
	    // if we're outside of 2 * map->radius, that should be safe
	    continue;
	  }
	  
	  for (int jz = -2; jz <= 2; jz++) {
	    for (int jy = -2; jy <= 2; jy++) {
	      for (int jx = -2; jx <= 2; jx++) {
		double xxx = MRC_MCRDX(crds, ix+jx, p);
		double yyy = MRC_MCRDY(crds, iy+jy, p);
		double zzz = MRC_MCRDZ(crds, iz+jz, p);
		double rrr = sqrt(sqr(xxx) + sqr(yyy) + sqr(zzz));
		if (rrr < r1) {
		  r1 -= map->dr;
		  goto loop2;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  r1 -= map->extra_dr * dr;

  MPI_Allreduce(&r1, &map->r1, 1, MPI_DOUBLE, MPI_MIN, ggcm_mhd_comm(mhd));
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup

void
ggcm_mhd_bnd_sphere_map_setup(struct ggcm_mhd_bnd_sphere_map *map, struct ggcm_mhd *mhd,
			      double radius, double dr, double extra_dr)
{
  map->mhd = mhd;
  map->radius = radius;
  map->dr = dr;
  map->extra_dr = extra_dr;
  ggcm_mhd_bnd_sphere_map_find_dr(map, &map->min_dr);
  ggcm_mhd_bnd_sphere_map_find_r1(map);

  ggcm_mhd_bnd_sphere_map_set_ymask(map);

  if (mhd->bnd_mask) {
    ggcm_mhd_bnd_sphere_map_set_bnd_mask(map);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_cc_n_map

static void
ggcm_mhd_bnd_sphere_map_find_cc_n_map(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  double r1 = map->r1;
  assert(r1 > 0.);

  int cc_n_map = 0;
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
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
	  double xx[3] = { MRC_MCRDX(crds, jx, p),
			   MRC_MCRDY(crds, jy, p),
			   MRC_MCRDZ(crds, jz, p), };
	  if (ggcm_mhd_bnd_sphere_map_is_ghost(map, xx)) {
	    cc_n_map++;
	  }
	}
      }
    }
  }
  map->cc_n_map = cc_n_map;
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_ec_n_map

static void
ggcm_mhd_bnd_sphere_map_find_ec_n_map(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  int ec_n_map[3] = {};
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);

    // edge-centered
    int l[3] = { 2, 2, 2 }, r[3] = { 1, 1, 1 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	l[d] = 0;
	r[d] = 0;
      }
    }
    for (int jz = -l[2]; jz < info.ldims[2] + r[2]; jz++) {
      for (int jy = -l[1]; jy < info.ldims[1] + r[1]; jy++) {
	for (int jx = -l[0]; jx < info.ldims[0] + r[0]; jx++) {
	  for (int d = 0; d < 3; d++) {
	    // find the correct edge centered coords for the locations of E
	    float crd_ec[3];
	    mrc_crds_at_ec(crds, jx,jy,jz, p, d, crd_ec);
	    double xx[3] = { crd_ec[0], crd_ec[1], crd_ec[2] };

	    if (ggcm_mhd_bnd_sphere_map_is_ghost(map, xx)) {
	      ec_n_map[d]++;
	    }
	  }
	}
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    map->ec_n_map[d] = ec_n_map[d];
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_fc_n_map
//
// NOTE: This is different from what we do for cc and ec:
// Here, we find the faces that make up the boundary between cells
// inside and outside the domain -- on these faces, we'll set
// reconstructed fluxes on the ghost side

static void
ggcm_mhd_bnd_sphere_map_find_fc_n_map(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  int fc_n_map[3] = {};
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);

    int l[3] = { 2, 2, 2 }, r[3] = { 1, 1, 1 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	l[d] = 0;
	r[d] = 0;
      }
    }
    for (int jz = -l[2]; jz < info.ldims[2] + r[2]; jz++) {
      for (int jy = -l[1]; jy < info.ldims[1] + r[1]; jy++) {
	for (int jx = -l[0]; jx < info.ldims[0] + r[0]; jx++) {
	  for (int d = 0; d < 3; d++) {
	    double xx0[3] = { MRC_MCRDX(crds, jx, p),
			      MRC_MCRDY(crds, jy, p),
			      MRC_MCRDZ(crds, jz, p) };
	    double xxp[3] = { MRC_MCRDX(crds, jx + (d == 0), p),
			      MRC_MCRDY(crds, jy + (d == 1), p),
			      MRC_MCRDZ(crds, jz + (d == 2), p) };
	    bool bnd0 = ggcm_mhd_bnd_sphere_map_is_bnd(map, xx0);
	    bool bndp = ggcm_mhd_bnd_sphere_map_is_bnd(map, xxp);
	    if (bnd0 != bndp) {
	      // one inside, the other outside
	      fc_n_map[d]++;
	    }
	  }
	}
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    map->fc_n_map[d] = fc_n_map[d];
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup_flds

void
ggcm_mhd_bnd_sphere_map_setup_flds(struct ggcm_mhd_bnd_sphere_map *map)
{
  // cell-centered
  ggcm_mhd_bnd_sphere_map_find_cc_n_map(map);
  mrc_fld_set_type(map->cc_imap, "int");
  mrc_fld_set_param_int_array(map->cc_imap, "dims", 2, (int[2]) { 4, map->cc_n_map });

  // edge-centered
  ggcm_mhd_bnd_sphere_map_find_ec_n_map(map);
  for (int d = 0; d < 3; d++) {
    mrc_fld_set_type(map->ec_imap[d], "int");
    mrc_fld_set_param_int_array(map->ec_imap[d], "dims", 2, (int[2]) { 4, map->ec_n_map[d] });
  }

  // face-centered
  ggcm_mhd_bnd_sphere_map_find_fc_n_map(map);
  for (int d = 0; d < 3; d++) {
    mrc_fld_set_type(map->fc_imap[d], "int");
    mrc_fld_set_param_int_array(map->fc_imap[d], "dims", 2, (int[2]) { 5, map->fc_n_map[d] });
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup_cc

void
ggcm_mhd_bnd_sphere_map_setup_cc(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  int cc_n_map = 0;
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
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
	  double xx[3] = { MRC_MCRDX(crds, jx, p),
			   MRC_MCRDY(crds, jy, p),
			   MRC_MCRDZ(crds, jz, p), };
	  if (ggcm_mhd_bnd_sphere_map_is_ghost(map, xx)) {
	    MRC_I2(map->cc_imap, 0, cc_n_map) = jx;
	    MRC_I2(map->cc_imap, 1, cc_n_map) = jy;
	    MRC_I2(map->cc_imap, 2, cc_n_map) = jz;
	    MRC_I2(map->cc_imap, 3, cc_n_map) = p;
	    cc_n_map++;
	  }
	}
      }
    }
  }

  assert(cc_n_map == map->cc_n_map);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup_ec

void
ggcm_mhd_bnd_sphere_map_setup_ec(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  int ec_n_map[3] = {};
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);

    int l[3] = { 2, 2, 2 }, r[3] = { 1, 1, 1 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	l[d] = 0;
	r[d] = 0;
      }
    }
    for (int jz = -l[2]; jz < info.ldims[2] + r[2]; jz++) {
      for (int jy = -l[1]; jy < info.ldims[1] + r[1]; jy++) {
	for (int jx = -l[0]; jx < info.ldims[0] + r[0]; jx++) {
	  for (int d = 0; d < 3; d++) {
	    // find the correct edge centered coords for the locations of E
	    float crd_ec[3];
	    mrc_crds_at_ec(crds, jx,jy,jz, p, d, crd_ec);
	    double xx[3] = { crd_ec[0], crd_ec[1], crd_ec[2] };

	    if (ggcm_mhd_bnd_sphere_map_is_ghost(map, xx)) {
	      MRC_I2(map->ec_imap[d], 0, ec_n_map[d]) = jx;
	      MRC_I2(map->ec_imap[d], 1, ec_n_map[d]) = jy;
	      MRC_I2(map->ec_imap[d], 2, ec_n_map[d]) = jz;
	      MRC_I2(map->ec_imap[d], 3, ec_n_map[d]) = p;
	      ec_n_map[d]++;
	    }
	  }
	}
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    assert(map->ec_n_map[d] == ec_n_map[d]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup_fc
//
// NOTE: This is different from what we do for cc and ec:
// Here, we find the faces that make up the boundary between cells
// inside and outside the domain -- on these faces, we'll set
// reconstructed fluxes on the ghost side

void
ggcm_mhd_bnd_sphere_map_setup_fc(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;

  struct mrc_fld *bnd_mask = NULL;
  if (mhd->bnd_mask) {
    bnd_mask = mrc_fld_get_as(mhd->bnd_mask, FLD_TYPE);
  }

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  int fc_n_map[3] = {};
  for (int p = 0; p < mrc_fld_nr_patches(mhd->fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);

    int l[3] = { 2, 2, 2 }, r[3] = { 1, 1, 1 };
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	l[d] = 0;
	r[d] = 0;
      }
    }
    for (int jz = -l[2]; jz < info.ldims[2] + r[2]; jz++) {
      for (int jy = -l[1]; jy < info.ldims[1] + r[1]; jy++) {
	for (int jx = -l[0]; jx < info.ldims[0] + r[0]; jx++) {
	  for (int d = 0; d < 3; d++) {
	    double xx0[3] = { MRC_MCRDX(crds, jx, p),
			      MRC_MCRDY(crds, jy, p),
			      MRC_MCRDZ(crds, jz, p) };
	    double xxp[3] = { MRC_MCRDX(crds, jx + (d == 0), p),
			      MRC_MCRDY(crds, jy + (d == 1), p),
			      MRC_MCRDZ(crds, jz + (d == 2), p) };
	    bool bnd0 = ggcm_mhd_bnd_sphere_map_is_bnd(map, xx0);
	    bool bndp = ggcm_mhd_bnd_sphere_map_is_bnd(map, xxp);
	    if (bnd0 != bndp) {
	      // one inside, the other outside
	      MRC_I2(map->fc_imap[d], 0, fc_n_map[d]) = jx + (d == 0);
	      MRC_I2(map->fc_imap[d], 1, fc_n_map[d]) = jy + (d == 1);
	      MRC_I2(map->fc_imap[d], 2, fc_n_map[d]) = jz + (d == 2);
	      MRC_I2(map->fc_imap[d], 3, fc_n_map[d]) = p;
	      MRC_I2(map->fc_imap[d], 4, fc_n_map[d]) = bndp;
	      fc_n_map[d]++;

	      if (bnd_mask) {
		// mark the interior cells adjacent to the boundary faces
		if (d == 0) {
		  M3(bnd_mask, 0, jx+bnd0,jy,jz, p) = 2.;
		} else if (d == 1) {
		  M3(bnd_mask, 0, jx,jy+bnd0,jz, p) = 2.;
		} else if (d == 2) {
		  M3(bnd_mask, 0, jx,jy,jz+bnd0, p) = 2.;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  for (int d = 0; d < 3; d++) {
    assert(map->fc_n_map[d] == fc_n_map[d]);
  }

  if (bnd_mask) {
    mrc_fld_put_as(bnd_mask, mhd->bnd_mask);
  }
}

