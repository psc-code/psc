
#include "mrc_ddc.h"

#include <mrc_domain.h>
#include <assert.h>
#include <stdio.h>

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
				   int dir[3], int off[3], int *p_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_local_patch_info(domain, p, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << pi.level, 1 << pi.level, 1 << pi.level };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = pi.idx3[d] + dir[d];
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

// this function incorporates Fujimoto (2011)-specific FDTD policy to
// some extent, ie., points on the boundary between coarse and fine levels
// are considered interior points on the coarse level, ghost points on the
// fine level.

bool
mrc_domain_is_ghost(struct mrc_domain *domain, int ext[3], int p, int i[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);

  int dir[3], dirx[3] = {};
  for (int d = 0; d < 3; d++) {
    if (i[d] < 0) {
      return true;
    } else if (ext[d] == 1 && i[d] == 0) {
      dir[d] = 0;
      dirx[d] = 1;
    } else if (i[d] < ldims[d]) {
      dir[d] = 0;
    } else if (ext[d] == 1 && i[d] == ldims[d]) {
      dir[d] = 1;
      dirx[d] = 1;
    } else {
      return true;
    }
  }
  // if outside, we've already returned true

  // inside, not on the boundary
  if (dir[0] == 0 && dirx[0] == 0 &&
      dir[1] == 0 && dirx[1] == 0) {
    return false;
  }

  // on the boundary...
  int dd[3];
  // do we border a coarse domain? (then it's a ghost point)
  for (dd[2] = 0; dd[2] >= 0; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	int p_nei;
	mrc_domain_get_neighbor_patch_coarse(domain, p, dd, &p_nei);
	if (p_nei >= 0) {
	  return true;
	}
      }
    }
  }

  // is another same level patch in line before us, then it's his, and we have
  // a ghost point
  for (dd[2] = 0; dd[2] >= 0; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	int p_nei;
	mrc_domain_get_neighbor_patch_same(domain, p, dd, &p_nei);
	if (p_nei >= 0) {
	  return p != p_nei;
	}
      }
    }
  }
  return true;
}

static void
mrc_domain_find_valid_point_same(struct mrc_domain *domain, int ext[3], int p, int i[3],
				 int *p_nei, int j[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);

  int dir[3], dirx[3] = {};
  for (int d = 0; d < 3; d++) {
    if (i[d] < 0) {
      dir[d] = -1;
    } else if (ext[d] == 1 && i[d] == 0) {
      dir[d] = 0;
      dirx[d] = 1;
    } else if (i[d] < ldims[d]) {
      dir[d] = 0;
    } else if (ext[d] == 1 && i[d] == ldims[d]) {
      dir[d] = 1;
      dirx[d] = 1;
    } else {
      dir[d] = 1;
    }
  }

  int dd[3];
  for (dd[2] = 0; dd[2] >= 0; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	mrc_domain_get_neighbor_patch_same(domain, p, dd, p_nei);
	if (*p_nei >= 0) {
	  for (int d = 0; d < 3; d++) {
	    j[d] = i[d] - dd[d] * ldims[d];
	  }
	  // need to double check whether we actually picked an interior point
	  if (!mrc_domain_is_ghost(domain, ext, *p_nei, j)) {
	    return;
	  }
	}
      }
    }
  }
  *p_nei = -1;
}

static void
mrc_domain_to_valid_point_same(struct mrc_domain *domain, int ext[3], int p, int i[3],
			       int *p_nei, int j[3])
{
  if (!mrc_domain_is_ghost(domain, ext, p, i)) {
    for (int d = 0; d < 3; d++) {
      j[d] = i[d];
    }
    *p_nei = p;
    return;
  }

  mrc_domain_find_valid_point_same(domain, ext, p, i, p_nei, j);
}

static void
mrc_domain_find_valid_point_coarse(struct mrc_domain *domain, int ext[3],
				   int p, int i[3], int *p_nei, int j[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  struct mrc_patch_info pi;
  mrc_domain_get_local_patch_info(domain, p, &pi);
    
  int ii[3], dir[3], dirx[3] = {};
  for (int d = 0; d < 3; d++) {
    ii[d] = i[d] + ((pi.idx3[d] & 1) ? ldims[d] / 2 : 0);

    if (ii[d] < 0) {
      dir[d] = -1;
    } else if (ext[d] == 1 && ii[d] == 0) {
      dir[d] = 0;
      dirx[d] = 1;
    } else if (ii[d] < ldims[d]) {
      dir[d] = 0;
    } else if (ext[d] == 1 && ii[d] == ldims[d]) {
      dir[d] = 1;
      dirx[d] = 1;
    } else {
      dir[d] = 1;
    }
  }

  int dd[3];
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	mrc_domain_get_neighbor_patch_coarse(domain, p, dd, p_nei);
	if (*p_nei >= 0) {
	  for (int d = 0; d < 3; d++) {
	    j[d] = ii[d] - dd[d] * ldims[d];
	  }
	  return;
	}
      }
    }
  }
  *p_nei = -1;
}

static void
mrc_domain_find_valid_point_fine(struct mrc_domain *domain, int p, int i[3],
				 int *p_nei, int j[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  struct mrc_patch_info pi;
  mrc_domain_get_local_patch_info(domain, p, &pi);
    
  int ii[3], off[3], dir[3];
  for (int d = 0; d < 3; d++) {
    off[d] = (i[d] + ldims[d]) / ldims[d] - 1;
    ii[d] = (i[d] + ldims[d]) % ldims[d];
    if (ii[d] < 0) {
      dir[d] = -1;
    } else if (ii[d] < ldims[d]) {
      dir[d] = 0;
    } else {
      dir[d] = 1;
    }
  }

  mrc_domain_get_neighbor_patch_fine(domain, p, dir, off, p_nei);
  if (*p_nei < 0) {
    return;
  }

  for (int d = 0; d < 3; d++) {
    j[d] = ii[d] - dir[d] * ldims[d];
  }
}

static inline int
div_2(int i)
{
  // divide by 2, but always round down
  return (i + 10) / 2 - 5;
}

static bool
mrc_ddc_amr_stencil_coarse(struct mrc_ddc *ddc, int ext[3],
			   struct mrc_ddc_amr_stencil *stencil,
			   int m, int p, int i[3])
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int p_nei, j[3];
  mrc_domain_find_valid_point_coarse(domain, ext, p,
				     (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
				     &p_nei, j);
  if (p_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil_entry *s = stencil->s; s < stencil->s + stencil->nr_entries; s++) {
    int jd[3], p_dnei, j_dnei[3];
    for (int d = 0; d < 3; d++) {
      jd[d] = j[d] + s->dx[d] * (i[d] & 1 && d < 2); // FIXME 3D
    }
    mrc_domain_to_valid_point_same(domain, ext, p_nei, jd, &p_dnei, j_dnei);
    assert(!mrc_domain_is_ghost(domain, ext, p_dnei, j_dnei));
    mrc_ddc_amr_add_value(ddc, p, m, i, p_dnei, m, j_dnei, s->val);
  }
  return true;
}

static bool
mrc_ddc_amr_stencil_fine(struct mrc_ddc *ddc, int ext[3],
			 struct mrc_ddc_amr_stencil *stencil,
			 int m, int p, int i[3])
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int p_nei, j[3];
  mrc_domain_find_valid_point_fine(domain, p, (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &p_nei, j);
  if (p_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil_entry *s = stencil->s; s < stencil->s + stencil->nr_entries; s++) {
    int id[3];
    for (int d = 0; d < 3; d++) {
      id[d] = 2*i[d] + s->dx[d];
    }
    mrc_domain_find_valid_point_fine(domain, p, id, &p_nei, j);
    assert(!mrc_domain_is_ghost(domain, ext, p_nei, j));
    mrc_ddc_amr_add_value(ddc, p, m, i, p_nei, m, j, s->val);
  }
  return true;
}

// ================================================================================

void
mrc_ddc_amr_set_by_stencil(struct mrc_ddc *ddc, int m, int bnd, int ext[3],
			   struct mrc_ddc_amr_stencil *stencil_coarse,
			   struct mrc_ddc_amr_stencil *stencil_fine)
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    int i[3];
    for (i[2] = 0; i[2] < ldims[2] + 0; i[2]++) { // FIXME 3D
      for (i[1] = -bnd; i[1] < ldims[1] + ext[1] + bnd; i[1]++) {
	for (i[0] = -bnd; i[0] < ldims[0] + ext[0] + bnd; i[0]++) {
	  if (i[0] >= ext[0] && i[0] < ldims[0] &&
	      i[1] >= ext[1] && i[1] < ldims[1] &&
	      i[2] >= ext[2] && i[2] < ldims[2]) {
	    assert(!mrc_domain_is_ghost(domain, ext, p, i));
	    continue;
	  }
	  if (!mrc_domain_is_ghost(domain, ext, p, i)) {
	    continue;
	  }

	  // at this point, we skipped all interior points, so only ghostpoints are left

	  // try to find an interior point corresponding to the current ghostpoint
	  int j[3], p_nei;
	  mrc_domain_find_valid_point_same(domain, ext, p, i, &p_nei, j);
	  if (p_nei >= 0) {
	    assert(!mrc_domain_is_ghost(domain, ext, p_nei, j));
	    mrc_ddc_amr_add_value(ddc, p, m, i, p_nei, m, j, 1.f);
	    continue;
	  }

	  // try to interpolate from coarse
	  if (mrc_ddc_amr_stencil_coarse(ddc, ext, stencil_coarse, m, p, i)) {
	    continue;
	  }
	      
	  // try to restrict from fine
	  if (mrc_ddc_amr_stencil_fine(ddc, ext, stencil_fine, m, p, i)) {
	    continue;
	  }

	  // oops, no way to fill this point?
	  // (This may be okay if the domain has holes or other physical boundaries)
	  MHERE;
	}
      }
    }
  }
}

