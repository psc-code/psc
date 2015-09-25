
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_bits.h>
#include <math.h>

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
// ggcm_mhd_bnd_sphere_map_find_r1_r2
//
// determines r1 to be small enough so that
// +/- 2 grid points around each cell with center outside of r2
// are outside (ie., their centers) the smaller r1 sphere

static void
ggcm_mhd_bnd_sphere_map_find_r1_r2(struct ggcm_mhd_bnd_sphere_map *map,
				   double radius, double *p_r1, double *p_r2)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double dr = map->min_dr;
  double r2 = radius;
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

  MPI_Allreduce(&r1, p_r1, 1, MPI_DOUBLE, MPI_MIN, ggcm_mhd_comm(mhd));

  *p_r2 = r2;
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_setup

void
ggcm_mhd_bnd_sphere_map_setup(struct ggcm_mhd_bnd_sphere_map *map, struct ggcm_mhd *mhd,
			      double radius)
{
  map->mhd = mhd;
  ggcm_mhd_bnd_sphere_map_find_dr(map, &map->min_dr);
  ggcm_mhd_bnd_sphere_map_find_r1_r2(map, radius, &map->r1, &map->r2);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_sphere_map_find_cc_n_map

void
ggcm_mhd_bnd_sphere_map_find_cc_n_map(struct ggcm_mhd_bnd_sphere_map *map)
{
  struct ggcm_mhd *mhd = map->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  double r1 = map->r1, r2 = map->r2;
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
  map->cc_n_map = cc_n_map;
}

