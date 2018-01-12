
#include "mrc_ddc.h"

#include <mrc_domain_private.h>
#include <assert.h>
#include <stdio.h>

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_same

void
mrc_domain_get_neighbor_patch_same(struct mrc_domain *domain, int gp,
				   int dx[3], int *gp_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_global_patch_info(domain, gp, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << pi.level, 1 << pi.level, 1 << pi.level };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = pi.idx3[d] + dx[d];
    if (domain->bc[d] == BC_PERIODIC && idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (domain->bc[d] == BC_PERIODIC && idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level, idx3, &pi_nei);
  *gp_nei = pi_nei.global_patch;
}

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_coarse

void
mrc_domain_get_neighbor_patch_coarse(struct mrc_domain *domain, int gp,
				     int dx[3], int *gp_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_global_patch_info(domain, gp, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << (pi.level - 1), 1 << (pi.level - 1), 1 << (pi.level - 1) };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = (pi.idx3[d] + dx[d] + 2) / 2 - 1;
    if (0&&idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (0&&idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level - 1, idx3, &pi_nei);
  *gp_nei = pi_nei.global_patch;
}

// ----------------------------------------------------------------------
// mrc_domain_get_neighbor_patch_fine

void
mrc_domain_get_neighbor_patch_fine(struct mrc_domain *domain, int gp,
				   int dir[3], int off[3], int *gp_nei)
{
  struct mrc_patch_info pi, pi_nei;
  mrc_domain_get_global_patch_info(domain, gp, &pi);
  // FIXME: how about if we only refine in selected directions?
  int mx[3] = { 1 << pi.level, 1 << pi.level, 1 << pi.level };
  int idx3[3];
  for (int d = 0; d < 3; d++) {
    idx3[d] = pi.idx3[d] + dir[d];
    if (0&&idx3[d] < 0) {
      idx3[d] += mx[d];
    }
    if (0&&idx3[d] >= mx[d]) {
      idx3[d] -= mx[d];
    }
    idx3[d] = 2 * idx3[d] + off[d];
  }
  mrc_domain_get_level_idx3_patch_info(domain, pi.level + 1, idx3, &pi_nei);
  *gp_nei = pi_nei.global_patch;
}

// ======================================================================

// this function incorporates Fujimoto (2011)-specific FDTD policy to
// some extent, ie., points on the boundary between coarse and fine levels
// are considered interior points on the coarse level, ghost points on the
// fine level.

bool
mrc_domain_is_ghost(struct mrc_domain *domain, int ext[3], int gp, int i[3])
{
#if 0 // FIXME, FDTD AMR relies on this definition of what's a ghost point
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);

  // FIXME simplify: dirx[d] == ext[d]
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
	int gp_nei;
	mrc_domain_get_neighbor_patch_coarse(domain, gp, dd, &gp_nei);
	if (gp_nei >= 0) {
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
	int gp_nei;
	mrc_domain_get_neighbor_patch_same(domain, gp, dd, &gp_nei);
	if (gp_nei >= 0) {
	  return gp != gp_nei;
	}
      }
    }
  }
  return true;
#else

  // flux correction
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
  //  mprintf("dir %d:%d dirx %d:%d\n", dir[0], dir[1], dirx[0], dirx[1]);
  // if outside, we've already returned true

  // inside, not on the boundary
  if (dir[0] == 0 && dirx[0] == 0 &&
      dir[1] == 0 && dirx[1] == 0 &&
      dir[2] == 0 && dirx[2] == 0) {
    return false;
  }

  // on the boundary...
  int dd[3];
  // do we border a fine domain? (then it's a ghost point)
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	int gp_nei;
	mrc_domain_get_neighbor_patch_fine(domain, gp, dd, (int [3]) { 0, 0, 0 }, &gp_nei);
	if (gp_nei >= 0) {
	  return true;
	}
      }
    }
  }

  // is another same level patch in line before us, then it's its, and we have
  // a ghost point
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	int gp_nei;
	mrc_domain_get_neighbor_patch_same(domain, gp, dd, &gp_nei);
	if (gp_nei >= 0) {
	  return gp != gp_nei;
	}
      }
    }
  }

  return false;
#endif
}

bool
mrc_domain_is_local_ghost(struct mrc_domain *domain, int ext[3], int lp, int i[3])
{
  struct mrc_patch_info pi;
  mrc_domain_get_local_patch_info(domain, lp, &pi);
  return mrc_domain_is_ghost(domain, ext, pi.global_patch,i);
}

void 
mrc_domain_find_valid_point_same(struct mrc_domain *domain, int ext[3], int gp, int i[3],
				 int *gp_nei, int j[3])
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
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  //	  continue;
	}
	mrc_domain_get_neighbor_patch_same(domain, gp, dd, gp_nei);
	if (*gp_nei >= 0) {
	  for (int d = 0; d < 3; d++) {
	    j[d] = i[d] - dd[d] * ldims[d];
	  }
	  // need to double check whether we actually picked an interior point
	  // FIXME!!! needs to use consistent _is_ghost()
	  if (!mrc_domain_is_ghost(domain, ext, *gp_nei, j)) {
	    return;
	  }
	}
      }
    }
  }
  *gp_nei = -1;
}

static void
mrc_domain_to_valid_point_same(struct mrc_domain *domain, int ext[3], int gp, int i[3],
			       int *gp_nei, int j[3])
{
  if (!mrc_domain_is_ghost(domain, ext, gp, i)) {
    for (int d = 0; d < 3; d++) {
      j[d] = i[d];
    }
    *gp_nei = gp;
    return;
  }

  mrc_domain_find_valid_point_same(domain, ext, gp, i, gp_nei, j);
}

void
mrc_domain_find_valid_point_coarse(struct mrc_domain *domain, int ext[3],
				   int gp, int i[3], int *gp_nei, int j[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  struct mrc_patch_info pi;
  mrc_domain_get_global_patch_info(domain, gp, &pi);
    
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
	mrc_domain_get_neighbor_patch_coarse(domain, gp, dd, gp_nei);
	if (*gp_nei >= 0) {
	  for (int d = 0; d < 3; d++) {
	    j[d] = ii[d] - dd[d] * ldims[d];
	  }
	  return;
	}
      }
    }
  }
  *gp_nei = -1;
}

void
mrc_domain_find_valid_point_fine(struct mrc_domain *domain, int ext[3], int gp, int i[3],
				 int *gp_nei, int j[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  struct mrc_patch_info pi;
  mrc_domain_get_global_patch_info(domain, gp, &pi);
    
  int off[3], dirl[3], dirh[3];
  for (int d = 0; d < 3; d++) {
    if (i[d] < 0) {
      dirl[d] = -1; dirh[d] = -1;
    } else if (i[d] == 0) {
      if (ext[d] == 1) { // on boundary
	dirl[d] = -1; dirh[d] = 0;
      } else {
	dirl[d] = 0; dirh[d] = 0;
      }
    } else if (i[d] < 2 * ldims[d]) {
        dirl[d] = 0; dirh[d] = 0;
    } else if (i[d] == 2 * ldims[d]) {
      if (ext[d] == 1) { // on boundary
	dirl[d] = 0; dirh[d] = 1;
      } else {
	dirl[d] = 1; dirh[d] = 1;
      }
    } else { // > 2 * ldims[d]
      dirl[d] = 1; dirh[d] = 1;
    }
  }

  bool found_somewhere = false;
  int dir[3];
  for (dir[2] = dirl[2]; dir[2] <= dirh[2]; dir[2]++) {
    for (dir[1] = dirl[1]; dir[1] <= dirh[1]; dir[1]++) {
      for (dir[0] = dirl[0]; dir[0] <= dirh[0]; dir[0]++) {
	for (int d = 0; d < 3; d++) {
	  if (dir[d] == -1) {
	    off[d] = 1;
	  } else if (dir[d] == 0) {
	    off[d] = (i[d] >= ldims[d]) ? 1 : 0;
	  } else { // dir[d] == 1
	    off[d] = 0;
	  }
	  j[d] = i[d] - 2 * ldims[d] * dir[d] - off[d] * ldims[d];
	}
	
	mrc_domain_get_neighbor_patch_fine(domain, gp, dir, off, gp_nei);
	if (*gp_nei >= 0) {
	  if (mrc_domain_is_ghost(domain, ext, *gp_nei, j)) {
	    found_somewhere = true;
	    continue;
	  }
	  return;
	}
      }
    }
  }
  assert(!found_somewhere);
  *gp_nei = -1;
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
			   int m, int gp, int i[3])
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int gp_nei, j[3];
  mrc_domain_find_valid_point_coarse(domain, ext, gp,
				     (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
				     &gp_nei, j);
  if (gp_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil_entry *s = stencil->s; s < stencil->s + stencil->nr_entries; s++) {
    int jd[3], gp_dnei, j_dnei[3];
    for (int d = 0; d < 3; d++) {
      jd[d] = j[d] + s->dx[d] * (i[d] & 1 && d < 2); // FIXME 3D
    }
    mrc_domain_to_valid_point_same(domain, ext, gp_nei, jd, &gp_dnei, j_dnei);
    assert(!mrc_domain_is_ghost(domain, ext, gp_dnei, j_dnei));
    mrc_ddc_amr_add_value(ddc, gp, m, i, gp_dnei, m, j_dnei, s->val);
  }
  return true;
}

static bool
mrc_ddc_amr_stencil_fine(struct mrc_ddc *ddc, int ext[3],
			 struct mrc_ddc_amr_stencil *stencil,
			 int m, int gp, int i[3])
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int gp_nei, j[3];
  struct mrc_patch_info info;
  mrc_domain_get_global_patch_info(domain, gp, &info);
  mrc_domain_find_valid_point_fine(domain, ext, gp, (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
  /* if (info.off[0] == 24 && gp_nei >= 0) { */
  /*   mprintf("gp %d off %d:%d:%d i %d:%d:%d gp_nei %d\n", gp, info.off[0], info.off[1], info.off[2], i[0], i[1], i[2], gp_nei); */
  /*   mrc_domain_get_global_patch_info(domain, gp_nei, &info); */
  /*   mprintf("off %d:%d:%d\n", info.off[0], info.off[1], info.off[2]); */
  /* } */
  if (gp_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil_entry *s = stencil->s; s < stencil->s + stencil->nr_entries; s++) {
    int id[3];
    for (int d = 0; d < 3; d++) {
      id[d] = 2*i[d] + s->dx[d];
    }
    mrc_domain_find_valid_point_fine(domain, ext, gp, id, &gp_nei, j);

    assert(!mrc_domain_is_ghost(domain, ext, gp_nei, j));
    mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, s->val);
  }
  return true;
}

// ================================================================================

void
mrc_ddc_amr_add_diagonal_one(struct mrc_ddc *ddc, int gp, int m, int i[3])
{
  mrc_ddc_amr_add_value(ddc, gp, m, i, gp, m, i, 1.f);
}

void
mrc_ddc_amr_set_by_stencil(struct mrc_ddc *ddc, int m, int bnd, int _ext[3],
			   struct mrc_ddc_amr_stencil *stencil_coarse,
			   struct mrc_ddc_amr_stencil *stencil_fine)
{
  struct mrc_domain *domain = mrc_ddc_get_domain(ddc);

  int ldims[3], gdims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  mrc_domain_get_global_dims(domain, gdims);
  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);

  int sw[3], ext[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = (_ext[d] == 0) ? bnd : 0;
    ext[d] = _ext[d];
    if (gdims[d] == 1) {
      sw[d] = 0;
      ext[d] = 0;
    }
  }
  for (int lp = 0; lp < nr_patches; lp++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(domain, lp, &info);
    int gp = info.global_patch;
    int i[3];
    for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
      for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	  if (i[0] >= ext[0] && i[0] < ldims[0] &&
	      i[1] >= ext[1] && i[1] < ldims[1] &&
	      i[2] >= ext[2] && i[2] < ldims[2]) {
	    assert(!mrc_domain_is_ghost(domain, ext, gp, i));
	    mrc_ddc_amr_add_diagonal_one(ddc, gp, m, i);
	    continue;
	  }

	  if (!mrc_domain_is_ghost(domain, ext, gp, i)) {
	    mrc_ddc_amr_add_diagonal_one(ddc, gp, m, i);
	    continue;
	  }

	  // at this point, we skipped all interior points, so only ghostpoints are left

	  // try to find an interior point corresponding to the current ghostpoint
	  int j[3], gp_nei;
	  mrc_domain_find_valid_point_same(domain, ext, gp, i, &gp_nei, j);
	  if (gp_nei >= 0) {
	    assert(!mrc_domain_is_ghost(domain, ext, gp_nei, j));
	    /* mprintf("gp %d i %d:%d:%d gp_nei %d j %d:%d:%d\n", */
	    /* 	    gp, i[0], i[1], i[2], gp_nei, j[0], j[1], j[2]); */
	    mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, 1.f);
	    continue;
	  }

	  // try to interpolate from coarse
	  if (stencil_coarse) {
	    if (mrc_ddc_amr_stencil_coarse(ddc, ext, stencil_coarse, m, gp, i)) {
	      continue;
	    }
	  }
	      
	  // try to restrict from fine
	  if (stencil_fine) {
	    if (mrc_ddc_amr_stencil_fine(ddc, ext, stencil_fine, m, gp, i)) {
	      continue;
	    }
	  }

	  // oops, no way to fill this point?
	  // (This may be okay if the domain has holes or other physical boundaries)
	  //MHERE;
	  mrc_ddc_amr_add_diagonal_one(ddc, gp, m, i);
	}
      }
    }
  }
}

