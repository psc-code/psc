
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_a3.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define AMR

// ----------------------------------------------------------------------
// Bz
// +-------+---+---+
// | o   o | x | x |
// |   X   +---O---+
// | o   o | x | x |
// +-------+---+---+

// Ey
// +-------+---+---+
// o   o   x   x   x
// X       X---+---O
// o   o   x   x   x
// +-------+---+---+

// Ez
// X---o---X---x---O
// |       |   |   |
// o   o   x---x---x
// |       |   |   |
// X---o---X---x---O

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_same

static void
mrc_domain_get_neighbor_patch_same(struct mrc_domain *domain, int p,
				   int dx[3], int *p_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_local_patch_info(domain, p, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << pi.level, 1 << pi.level, 1 << pi.level };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = pi.idx3[d] + dx[d];
    if (idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level, idx3, &pi_nei);
  *p_nei = pi_nei.patch;
}

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_coarse

static void
mrc_domain_get_neighbor_patch_coarse(struct mrc_domain *domain, int p,
				     int dx[3], int *p_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_local_patch_info(domain, p, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << (pi.level - 1), 1 << (pi.level - 1), 1 << (pi.level - 1) };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = (pi.idx3[d] + dx[d] + 2) / 2 - 1;
    if (idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level - 1, idx3, &pi_nei);
  *p_nei = pi_nei.patch;
}

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_fine

static void
mrc_domain_get_neighbor_patch_fine(struct mrc_domain *domain, int p,
				   int dx[3], int off[3], int *p_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_local_patch_info(domain, p, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << pi.level, 1 << pi.level, 1 << pi.level };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = pi.idx3[d] + dx[d];
    if (idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
    idx3[d] = 2 * idx3[d] + off[d];
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level + 1, idx3, &pi_nei);
  *p_nei = pi_nei.patch;
}

// ======================================================================

enum {
  EX,
  EY,
  EZ,
  HX,
  HY,
  HZ,
  NR_COMPS,
};

// ======================================================================

static void
copy_block(struct mrc_a3_patch *fldp_to, struct mrc_a3_patch *fldp_from,
	   int off_to[3], int off_from[3], int len[3], int m)
{
  for (int iz = 0; iz < len[2]; iz++) {
    for (int iy = 0; iy < len[1]; iy++) {
      for (int ix = 0; ix < len[0]; ix++) {
	MRC_A3(fldp_to, m, ix + off_to[0], iy + off_to[1], iz + off_to[2]) =
	  MRC_A3(fldp_from, m, ix + off_from[0], iy + off_from[1], iz + off_from[2]);
      }
    }
  }
}


static void
fill_ghosts_same(struct mrc_a3 *fld, int m, int bnd, int ext[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    int p_nei;
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  mrc_domain_get_neighbor_patch_same(fld->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
	    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
	    struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
	    int off_to[3], off_from[3], len[3];
	    for (int d = 0; d < 3; d++) {
	      if (dir[d] < 0) {
		off_to[d]   = -bnd;
		// FIXME we shouldn't extend the bnd region in general, definitely not into our own
		// points, really only if there's a fine grid which means there's no other way to
		// set the point(s)
		len[d]      = bnd + (dir[1-d] != 0 ? ext[d] : 0);
	      } else if (dir[d] == 0) {
		off_to[d]   = 0;
		len[d]      = ext[d] + ldims[d];
	      } else {
		off_to[d]   = ldims[d];
		len[d]      = ext[d] + bnd;
	      }
	      off_from[d] = off_to[d] - dir[d] * ldims[d];
	    }
	    copy_block(fldp, fldp_nei, off_to, off_from, len, m);
	  }
	}
      }
    }
  }
}

static void
fill_ghosts_H_same(struct mrc_a3 *fld)
{
  fill_ghosts_same(fld, HZ, 2, (int[]) { 0, 0, 1 });
}

static void
fill_ghosts_E_same(struct mrc_a3 *fld)
{
  fill_ghosts_same(fld, EY, 2, (int[]) { 1, 0, 1 });
  fill_ghosts_same(fld, EZ, 2, (int[]) { 1, 1, 0 });
}

static void
fill_ghosts_H_coarse(struct mrc_a3 *fld)
{
  int bnd = 2, ext[3] = { 0, 0, 1 };
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(fld->domain, p, &pi);
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int p_nei;
	  mrc_domain_get_neighbor_patch_coarse(fld->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
	    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
	    struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
	    int off_from[3], off_to[3], len[3];
	    for (int d = 0; d < 3; d++) {
	      if (dir[d] == -1) {
		len[d] = bnd;
		off_to[d] = -bnd;
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = -bnd + 2*ldims[d];
		} else {
		  off_from[d] = -bnd + ldims[d];
		}
	      } else if (dir[d] == 0) {
		len[d] = ldims[d];
		off_to[d] = 0;
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = 0;
		} else {
		  off_from[d] = ldims[d];
		}
	      } else {
		len[d] = bnd;
		off_to[d] = ext[d] + ldims[d];
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = ext[d] + ldims[d];
		} else {
		  off_from[d] = ext[d];
		}
	      }
	    }
	    for (int iz = 0; iz < len[2]; iz++) {
	      for (int iy = 0; iy < len[1]; iy++) {
		for (int ix = 0; ix < len[0]; ix++) {
		  MRC_A3(fldp, HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2]) =
		    MRC_A3(fldp_nei, HZ, (ix + off_from[0])/2, (iy + off_from[1])/2, (iz + off_from[2])/2);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

static void
fill_ghosts_E_coarse(struct mrc_a3 *fld)
{
  int bnd = 2, ext[3] = { 1, 0, 1 };
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(fld->domain, p, &pi);
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int p_nei;
	  mrc_domain_get_neighbor_patch_coarse(fld->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
	    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
	    struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
	    int off_from[3], off_to[3], len[3];
	    for (int d = 0; d < 3; d++) {
	      if (dir[d] == -1) {
		len[d] = bnd + ext[d];
		off_to[d] = -bnd;
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = -bnd + 2*ldims[d];
		} else {
		  off_from[d] = -bnd + ldims[d];
		}
	      } else if (dir[d] == 0) {
		len[d] = ldims[d];
		off_to[d] = 0;
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = 0;
		} else {
		  off_from[d] = ldims[d];
		}
		len[d] += ext[d];
	      } else {
		len[d] = bnd + ext[d];
		off_to[d] = ldims[d];
		if ((pi.idx3[d] & 1) == 0) {
		  off_from[d] = ldims[d];
		} else {
		  off_from[d] = 0;
		}
	      }
	    }
	    for (int iz = 0; iz < len[2]; iz++) {
	      for (int iy = 0; iy < len[1]; iy++) {
		for (int ix = 0; ix < len[0]; ix++) {
		  int k = (ix + off_to[0]) & 1;
		  MRC_A3(fldp, EY, ix + off_to[0], iy + off_to[1], iz + off_to[2]) =
		    .5f * (MRC_A3(fldp_nei, EY, (ix + off_from[0])/2    , (iy + off_from[1])/2, iz) + 
			   MRC_A3(fldp_nei, EY, (ix + off_from[0])/2 + k, (iy + off_from[1])/2, iz));
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

static void
fill_ghosts_H_fine(struct mrc_a3 *fld)
{
  int bnd = 2;
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);

    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int sub[3];
	  int len[3];
	  for (int d = 0; d < 3; d++) {
	    if (dir[d] == 0) {
	      sub[d] = 2;
	    } else {
	      sub[d] = 1;
	    }
	  }
	  int is[3];
	  for (is[2] = 0; is[2] < sub[2]; is[2]++) {
	    for (is[1] = 0; is[1] < sub[1]; is[1]++) {
	      for (is[0] = 0; is[0] < sub[0]; is[0]++) {
		int off[3], off_to[3], off_from[3];
		for (int d = 0; d < 3; d++) {
		  if (dir[d] == -1) {
		    off[d] = 1;
		    len[d] = bnd;
		    off_to[d] = -bnd;
		    off_from[d] = 2 * -bnd + ldims[d];
		  } else if (dir[d] == 0) {
		    off[d] = is[d];
		    len[d] = ldims[d] / 2;
		    off_to[d] = is[d] * len[d];
		    off_from[d] = 0;
		  } else {
		    off[d] = 0;
		    len[d] = bnd;
		    off_to[d] = ldims[d];
		    off_from[d] = 0;
		  }
		}

		int p_nei;
		mrc_domain_get_neighbor_patch_fine(fld->domain, p, dir, off, &p_nei);
		if (p_nei < 0) {
		  continue;
		}
		struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);

		for (int iz = 0; iz < 1; iz++) {
		  for (int iy = 0; iy < len[1]; iy++) {
		    for (int ix = 0; ix < len[0]; ix++) {
		      MRC_A3(fldp, HZ, ix + off_to[0],iy + off_to[1],iz) =
			(1./8.f) * (2.f * MRC_A3(fldp_nei, HZ, 2*ix   + off_from[0],2*iy   + off_from[1],iz) +
				    2.f * MRC_A3(fldp_nei, HZ, 2*ix+1 + off_from[0],2*iy   + off_from[1],iz) +
				    2.f * MRC_A3(fldp_nei, HZ, 2*ix   + off_from[0],2*iy+1 + off_from[1],iz) +
				    2.f * MRC_A3(fldp_nei, HZ, 2*ix+1 + off_from[0],2*iy+1 + off_from[1],iz));
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

static void
fill_ghosts_E_fine(struct mrc_a3 *fld)
{
  int bnd = 2, ext[3] = { 1, 0, 1 };
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);

    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int sub[3];
	  int len[3];
	  for (int d = 0; d < 3; d++) {
	    if (dir[d] == 0) {
	      sub[d] = 2;
	    } else {
	      sub[d] = 1;
	    }
	  }
	  int is[3];
	  for (is[2] = 0; is[2] < sub[2]; is[2]++) {
	    for (is[1] = 0; is[1] < sub[1]; is[1]++) {
	      for (is[0] = 0; is[0] < sub[0]; is[0]++) {
		int off[3], off_to[3], off_from[3];
		for (int d = 0; d < 3; d++) {
		  if (dir[d] == -1) {
		    off[d] = 1;
		    len[d] = bnd;
		    off_to[d] = -bnd;
		    off_from[d] = 2 * -bnd + ldims[d];
		  } else if (dir[d] == 0) {
		    off[d] = is[d];
		    len[d] = ldims[d] / 2 - (is[d] == 0 ? ext[d] : 0);
		    off_to[d] = is[d] * len[d] + (is[d] == 0 ? ext[d] : 0);
		    off_from[d] = 2 * (is[d] == 0 ? ext[d] : 0);
		  } else {
		    off[d] = 0;
		    len[d] = bnd;
		    off_to[d] = ext[d] + ldims[d];
		    off_from[d] = 2 * ext[d];
		  }
		}
		int p_nei;
		mrc_domain_get_neighbor_patch_fine(fld->domain, p, dir, off, &p_nei);
		if (p_nei < 0) {
		  continue;
		}

		struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
		for (int iz = 0; iz < 1; iz++) {
		  for (int iy = 0; iy < len[1]; iy++) {
		    for (int ix = 0; ix < len[0]; ix++) {
		      MRC_A3(fldp, EY, ix + off_to[0],iy + off_to[1],iz) =
			(p == 11 ? 0. : 0.) +
			(1./8.f) * (2.f * MRC_A3(fldp_nei, EY, 2*ix   + off_from[0],2*iy   + off_from[1],iz) +
				    2.f * MRC_A3(fldp_nei, EY, 2*ix   + off_from[0],2*iy+1 + off_from[1],iz) +
				    1.f * MRC_A3(fldp_nei, EY, 2*ix+1 + off_from[0],2*iy   + off_from[1],iz) +
				    1.f * MRC_A3(fldp_nei, EY, 2*ix+1 + off_from[0],2*iy+1 + off_from[1],iz) +
				    1.f * MRC_A3(fldp_nei, EY, 2*ix-1 + off_from[0],2*iy   + off_from[1],iz) +
				    1.f * MRC_A3(fldp_nei, EY, 2*ix-1 + off_from[0],2*iy+1 + off_from[1],iz));
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

#define F3 MRC_A3

// ======================================================================

static void
fill_ghosts_H(struct mrc_a3 *fld)
{
  fill_ghosts_H_same(fld);
  fill_ghosts_H_coarse(fld); // not needed for fdtd
  fill_ghosts_H_fine(fld);
}

static void
fill_ghosts_E(struct mrc_a3 *fld)
{
  fill_ghosts_E_same(fld);
  fill_ghosts_E_coarse(fld);
  fill_ghosts_E_fine(fld); // not needed for fdtd
}

static void
fill_ghosts(struct mrc_a3 *fld)
{
  fill_ghosts_H(fld);
  fill_ghosts_E(fld);
}

// ----------------------------------------------------------------------
// step_fdtd

static void
step_fdtd(struct mrc_a3 *fld)
{
  struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
#ifdef AMR
  float dt = 1. / 32;
#else
  float dt = 1. / 16;
#endif

  fill_ghosts_H(fld);

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, EX, ix,iy,iz) +=
      	cny * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix,iy-1,iz)) -
      	cnz * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix,iy,iz-1));

      F3(fldp, EY, ix,iy,iz) +=
	cnz * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy,iz-1)) -
	cnx * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix-1,iy,iz));

      F3(fldp, EZ, ix,iy,iz) +=
      	cnx * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix-1,iy,iz)) -
      	cny * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy-1,iz));
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  fill_ghosts_E(fld);

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 0) {
      F3(fldp, HX, ix,iy,iz) -=
	cny * (F3(fldp, EZ, ix,iy+1,iz) - F3(fldp, EZ, ix,iy,iz)) -
	cnz * (F3(fldp, EY, ix,iy,iz+1) - F3(fldp, EY, ix,iy,iz));
      
      F3(fldp, HY, ix,iy,iz) -=
	cnz * (F3(fldp, EX, ix,iy,iz+1) - F3(fldp, EX, ix,iy,iz)) -
	cnx * (F3(fldp, EZ, ix+1,iy,iz) - F3(fldp, EZ, ix,iy,iz));
      
      F3(fldp, HZ, ix,iy,iz) -=
	cnx * (F3(fldp, EY, ix+1,iy,iz) - F3(fldp, EY, ix,iy,iz)) -
	cny * (F3(fldp, EX, ix,iy+1,iz) - F3(fldp, EX, ix,iy,iz));
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 0) {
      F3(fldp, HX, ix,iy,iz) -=
	cny * (F3(fldp, EZ, ix,iy+1,iz) - F3(fldp, EZ, ix,iy,iz)) -
	cnz * (F3(fldp, EY, ix,iy,iz+1) - F3(fldp, EY, ix,iy,iz));
      
      F3(fldp, HY, ix,iy,iz) -=
	cnz * (F3(fldp, EX, ix,iy,iz+1) - F3(fldp, EX, ix,iy,iz)) -
	cnx * (F3(fldp, EZ, ix+1,iy,iz) - F3(fldp, EZ, ix,iy,iz));
      
      F3(fldp, HZ, ix,iy,iz) -=
	cnx * (F3(fldp, EY, ix+1,iy,iz) - F3(fldp, EY, ix,iy,iz)) -
	cny * (F3(fldp, EX, ix,iy+1,iz) - F3(fldp, EX, ix,iy,iz));
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  fill_ghosts_H(fld);

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, EX, ix,iy,iz) +=
	cny * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix,iy-1,iz)) -
	cnz * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix,iy,iz-1));
      
      F3(fldp, EY, ix,iy,iz) +=
	cnz * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy,iz-1)) -
	cnx * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix-1,iy,iz));
      
      F3(fldp, EZ, ix,iy,iz) +=
	cnx * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix-1,iy,iz)) -
	cny * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy-1,iz));
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  mrc_domain_set_type(domain, "amr");
  mrc_domain_set_param_int3(domain, "m", (int [3]) { 8, 8, 1});
  mrc_crds_set_type(crds, "amr_uniform");
  mrc_crds_set_param_int(crds, "sw", 3);
  
  mrc_domain_set_from_options(domain);
#if 0
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
#else
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
#ifdef AMR
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
#else
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
#endif
#endif

  mrc_domain_setup(domain);
  mrc_domain_plot(domain);

  // create and fill a field

  struct mrc_a3 *fld = mrc_domain_a3_create(domain);
  mrc_a3_set_name(fld, "fld");
  mrc_a3_set_param_int(fld, "nr_comps", NR_COMPS);
  mrc_a3_set_param_int(fld, "sw", 3);
  mrc_a3_set_from_options(fld);
  mrc_a3_setup(fld);
  mrc_a3_set_comp_name(fld, EX, "EX");
  mrc_a3_set_comp_name(fld, EY, "EY");
  mrc_a3_set_comp_name(fld, EZ, "EZ");
  mrc_a3_set_comp_name(fld, HX, "HX");
  mrc_a3_set_comp_name(fld, HY, "HY");
  mrc_a3_set_comp_name(fld, HZ, "HZ");

  float kx = 2. * M_PI, ky = 2. * M_PI;

  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);

#if 1
    mrc_a3_foreach(fldp, ix,iy,iz, 3, 3) {
      MRC_A3(fldp, EY, ix,iy,iz) = 1.f / 0.f;
      MRC_A3(fldp, HZ, ix,iy,iz) = 1.f / 0.f;
    } mrc_a3_foreach_end;
#endif
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 1) {
      float x_cc = MRC_MCRDX(crds, ix);
      float y_cc = MRC_MCRDY(crds, iy);
      float x_nc = .5f * (MRC_MCRDX(crds, ix-1) + MRC_MCRDX(crds, ix));
      float y_nc = .5f * (MRC_MCRDY(crds, iy-1) + MRC_MCRDY(crds, iy));
      //      MRC_A3(fldp, EZ, ix,iy,iz) = sin(.5+kx * x_nc) * cos(.5+ky * y_nc);
      if (iy < ldims[1]) {
	MRC_A3(fldp, EY, ix,iy,iz) = sin(.5+kx * x_nc);// * cos(.5+ky * y_cc);
	if (ix < ldims[0]) {
	  MRC_A3(fldp, HZ, ix,iy,iz) = sin(.5+kx * x_cc);// * cos(.5+ky * y_cc);
	}
      }
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  // write field to disk

  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_type(io, "ascii");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  for (int n = 0; n <= 32; n++) {
    fill_ghosts(fld); // FIXME, only for output/debug

    mrc_io_open(io, "w", n, n);
    mrc_a3_write(fld, io);
    mrc_io_close(io);

    step_fdtd(fld);
  }

  mrc_io_destroy(io);

  mrc_a3_destroy(fld);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
