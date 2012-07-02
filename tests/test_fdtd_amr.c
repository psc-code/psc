
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
// |     o | x | x |
// |   X   +---O---+
// |     o | x | x |
// +-------+---+---+

// Ey
// +-------+---+---+
// o   o   x   x   x
// X       X---+---O
// o   o   x   x   x
// +-------+---+---+

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
    if ((dx[d] == -1 && (pi.idx3[d] & 1) == 1) ||
	(dx[d] ==  1 && (pi.idx3[d] & 1) == 0)) {
      *p_nei = -1;
      return;
    }

    idx3[d] = pi.idx3[d] / 2 + dx[d];
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
fill_ghosts_H_same(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    int p_nei;
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { -1, 0, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	for (int iy = 0; iy < ldims[1]; iy++) {
	  MRC_A3(fldp, HZ, -1,iy,iz) = MRC_A3(fldp_nei, HZ, ldims[0]-1,iy,iz);
	}
      }
    }
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { 0, -1, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	for (int ix = 0; ix < ldims[0]; ix++) {
	  MRC_A3(fldp, HZ, ix,-1,iz) = MRC_A3(fldp_nei, HZ, ix,ldims[1]-1,iz);
	}
      }
    }
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { -1, -1, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	MRC_A3(fldp, HZ, -1,-1,iz) = MRC_A3(fldp_nei, HZ, ldims[0]-1,ldims[1]-1,iz);
      }
    }
  }
}

static void
fill_ghosts_E_same(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    int p_nei;
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { 1, 0, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	for (int iy = 0; iy < ldims[1]; iy++) {
	  MRC_A3(fldp, EY, ldims[0],iy,iz) = MRC_A3(fldp_nei, EY, 0,iy,iz);
	}
      }
    }
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { 0, 1, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	for (int ix = 0; ix < ldims[0]; ix++) {
	  MRC_A3(fldp, EY, ix,ldims[1],iz) = MRC_A3(fldp_nei, EY, ix,0,iz);
	}
      }
    }
    mrc_domain_get_neighbor_patch_same(fld->domain, p, (int [3]) { 1, 1, 0 }, &p_nei);
    if (p_nei >= 0) {
      struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
      struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
      for (int iz = 0; iz < ldims[2]; iz++) {
	MRC_A3(fldp, EY, ldims[0],ldims[1],iz) = MRC_A3(fldp_nei, EY, 0,0,iz);
      }
    }
  }
}

static void
fill_ghosts_H_coarse(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(fld->domain, p, &pi);
    int p_nei;
    int dx[3] = { -1, 0, 0 };
    mrc_domain_get_neighbor_patch_coarse(fld->domain, p, dx, &p_nei);
    if (p_nei < 0) {
      continue;
    }
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
    for (int iz = 0; iz < ldims[2]; iz++) {
      for (int iy = 0; iy < ldims[1]; iy++) {
	int iy_nei = (iy + ((pi.idx3[1] & 1) == 1 ? ldims[1] : 0)) / 2;
	MRC_A3(fldp, HZ, -1,iy,iz) = MRC_A3(fldp_nei, HZ, ldims[0]-1,iy_nei,iz);
      }
    }
  }
}

static void
fill_ghosts_E_coarse(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(fld->domain, p, &pi);
    int p_nei;
    int dx[3] = { 1, 0, 0 };
    mrc_domain_get_neighbor_patch_coarse(fld->domain, p, dx, &p_nei);
    if (p_nei < 0) {
      continue;
    }
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    struct mrc_a3_patch *fldp_nei = mrc_a3_patch_get(fld, p_nei);
    for (int iz = 0; iz < ldims[2]; iz++) {
      for (int iy = 0; iy < ldims[1]; iy++) {
	int iy_nei = (iy + ((pi.idx3[1] & 1) == 1 ? ldims[1] : 0)) / 2;
	int j = iy & 1;
	MRC_A3(fldp, EY, ldims[0],iy,iz) =
	  .5f * (MRC_A3(fldp_nei, EY, 0,iy_nei  ,iz) +
		 MRC_A3(fldp_nei, EY, 0,iy_nei+j,iz));
      }
    }
  }
}

static void
fill_ghosts_H_fine(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    int p_nei1, p_nei2;
    int dx[3] = { -1, 0, 0 };
    int off[3] = { 1, 0, 0 };
    mrc_domain_get_neighbor_patch_fine(fld->domain, p, dx, off, &p_nei1);
    if (p_nei1 < 0) {
      continue;
    }
    off[1] = 1;
    mrc_domain_get_neighbor_patch_fine(fld->domain, p, dx, off, &p_nei2);
    assert(p_nei2 >= 0);
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    struct mrc_a3_patch *fldp_nei1 = mrc_a3_patch_get(fld, p_nei1);
    struct mrc_a3_patch *fldp_nei2 = mrc_a3_patch_get(fld, p_nei2);
    for (int iz = 0; iz < ldims[2]; iz++) {
      for (int iy = 0; iy < ldims[1] / 2; iy++) {
	MRC_A3(fldp, HZ, -1,iy,iz) =
	  (1./8.f) * (2.f * MRC_A3(fldp_nei1, HZ, ldims[0]-2,2*iy  ,iz) +
		      2.f * MRC_A3(fldp_nei1, HZ, ldims[0]-1,2*iy  ,iz) +
		      2.f * MRC_A3(fldp_nei1, HZ, ldims[0]-2,2*iy+1,iz) +
		      2.f * MRC_A3(fldp_nei1, HZ, ldims[0]-1,2*iy+1,iz));
      }
      for (int iy = 0; iy < ldims[1] / 2; iy++) {
	MRC_A3(fldp, HZ, -1,iy + ldims[1] / 2,iz) =
	  (1./8.f) * (2.f * MRC_A3(fldp_nei2, HZ, ldims[0]-2,2*iy  ,iz) +
		      2.f * MRC_A3(fldp_nei2, HZ, ldims[0]-1,2*iy  ,iz) +
		      2.f * MRC_A3(fldp_nei2, HZ, ldims[0]-2,2*iy+1,iz) +
		      2.f * MRC_A3(fldp_nei2, HZ, ldims[0]-1,2*iy+1,iz));
      }
    }
  }
}

static void
fill_ghosts_E_fine(struct mrc_a3 *fld)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  for (int p = 0; p < fld->nr_patches; p++) {
    int p_nei1, p_nei2;
    int dx[3] = { 1, 0, 0 };
    int off[3] = { 0, 0, 0 };
    mrc_domain_get_neighbor_patch_fine(fld->domain, p, dx, off, &p_nei1);
    if (p_nei1 < 0) {
      continue;
    }
    off[1] = 1;
    mrc_domain_get_neighbor_patch_fine(fld->domain, p, dx, off, &p_nei2);
    assert(p_nei2 >= 0);
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    struct mrc_a3_patch *fldp_nei1 = mrc_a3_patch_get(fld, p_nei1);
    struct mrc_a3_patch *fldp_nei2 = mrc_a3_patch_get(fld, p_nei2);
#if 0
    for (int iz = 0; iz < ldims[2]; iz++) {
      for (int iy = 0; iy < ldims[1] / 2; iy++) {
	MRC_A3(fldp, EY, ldims[0],iy,iz) = 
	  (1.f/4.f) * (1.f  * MRC_A3(fldp_nei1, EY, 0,2*iy  ,iz) +
		       .5f  * MRC_A3(fldp_nei1, EY, 1,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei1, EY,-1,2*iy+1,iz) + // FIXME is -1 valid?!
		       1.f  * MRC_A3(fldp_nei1, EY, 0,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei1, EY, 1,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei1, EY,-1,2*iy+1,iz));
      }
      for (int iy = 0; iy < ldims[1] / 2; iy++) {
	MRC_A3(fldp, EY, ldims[0],iy + ldims[1] / 2,iz) =
	  (1.f/4.f) * (1.f  * MRC_A3(fldp_nei2, EY, 0,2*iy  ,iz) +
		       .5f  * MRC_A3(fldp_nei2, EY, 1,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei2, EY,-1,2*iy+1,iz) + // FIXME is -1 valid?!
		       1.f  * MRC_A3(fldp_nei2, EY, 0,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei2, EY, 1,2*iy+1,iz) +
		       .5f  * MRC_A3(fldp_nei2, EY,-1,2*iy+1,iz));
      }
    }
#endif
  }
}

#define F3 MRC_A3

// ======================================================================

static void
fill_ghosts_H(struct mrc_a3 *fld)
{
  fill_ghosts_H_same(fld);
  fill_ghosts_H_coarse(fld);
  fill_ghosts_H_fine(fld);
}

static void
fill_ghosts_E(struct mrc_a3 *fld)
{
  fill_ghosts_E_same(fld);
  fill_ghosts_E_coarse(fld);
  fill_ghosts_E_fine(fld);
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
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 0) {
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
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 0) {
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
  mrc_crds_set_param_int(crds, "sw", 1);
  
  mrc_domain_set_from_options(domain);
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

  mrc_domain_setup(domain);
  mrc_domain_plot(domain);

  // create and fill a field

  struct mrc_a3 *fld = mrc_domain_a3_create(domain);
  mrc_a3_set_name(fld, "fld");
  mrc_a3_set_param_int(fld, "nr_comps", NR_COMPS);
  mrc_a3_set_param_int(fld, "sw", 1);
  mrc_a3_set_from_options(fld);
  mrc_a3_setup(fld);
  mrc_a3_set_comp_name(fld, EX, "EX");
  mrc_a3_set_comp_name(fld, EY, "EY");
  mrc_a3_set_comp_name(fld, EZ, "EZ");
  mrc_a3_set_comp_name(fld, HX, "HX");
  mrc_a3_set_comp_name(fld, HY, "HY");
  mrc_a3_set_comp_name(fld, HZ, "HZ");

  float kx = 2. * M_PI, ky = 2. * M_PI;

  mrc_a3_foreach_patch(fld, p) {
    struct mrc_a3_patch *fldp = mrc_a3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    mrc_a3_foreach(fldp, ix,iy,iz, 0, 1) {
      float x_cc = MRC_MCRDX(crds, ix);
      float y_cc = MRC_MCRDY(crds, iy);
      float x_nc = .5f * (MRC_MCRDX(crds, ix-1) + MRC_MCRDX(crds, ix));
      //float y_nc = .5f * (MRC_MCRDY(crds, iy-1) + MRC_MCRDY(crds, iy));
      MRC_A3(fldp, EY, ix,iy,iz) = sin(.5+kx * x_nc) * cos(.5+ky * y_cc);
      MRC_A3(fldp, HZ, ix,iy,iz) = sin(.5+kx * x_cc) * cos(.5+ky * y_cc);
    } mrc_a3_foreach_end;
    mrc_a3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  // write field to disk

  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_type(io, "ascii");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  for (int n = 0; n < 10; n++) {
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
