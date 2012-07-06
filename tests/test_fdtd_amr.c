
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define AMR

// ----------------------------------------------------------------------

// Bx
// +-------+---+---+
// |   o   o   x   x
// X       X---+---O
// |   o   o   x   x
// +-------+---+---+

// By
// +---X-o-+-x-O-x-+
// |       |   |   |
// |     o +-x-+-x-+
// |       |   |   |
// +---X-o-+-x-O-x-+

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

// FIXME
const int max_rows = 10000;
const int max_entries = 20000;

struct mrc_ddc_amr_row {
  int patch;
  int idx;
  int first_entry;
};

struct mrc_ddc_amr_entry {
  int patch;
  int idx;
  float val;
};

struct mrc_ddc_amr {
  struct mrc_ddc_amr_row *rows;
  struct mrc_ddc_amr_entry *entries;

  int nr_rows;
  int nr_entries;

  struct mrc_domain *domain;
  int sw;
  int ib[3], im[3];
};

// ----------------------------------------------------------------------
// mrc_ddc_amr_setup

static void
mrc_ddc_amr_setup(struct mrc_ddc_amr *amr)
{
  assert(amr->domain);
  assert(amr->sw);

  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  // needs to be compatible with how mrc_m3 indexes
  for (int d = 0; d < 3; d++) {
    amr->ib[d] = -amr->sw;
    amr->im[d] = ldims[d] + 2 * amr->sw;
  }

  amr->rows = calloc(max_rows + 1, sizeof(*amr->rows));
  amr->entries = calloc(max_entries, sizeof(*amr->entries));

  amr->nr_entries = 0;
  amr->nr_rows = 0;
}

// ----------------------------------------------------------------------
// mrc_ddc_add_destroy

static void
mrc_ddc_amr_destroy(struct mrc_ddc_amr *amr)
{
  free(amr->rows);
  free(amr->entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_add_value

static void
mrc_ddc_amr_add_value(struct mrc_ddc_amr *amr,
		      int row_patch, int rowm, int rowx, int rowy, int rowz,
		      int col_patch, int colm, int colx, int coly, int colz,
		      float val)
{
  // FIXME, a * F3(i,j,k) + b * F3(i,j,k) should be combined as (a + b) * F3(i,j,k)
  // WARNING, all elements for any given row must be added contiguously!

  assert(row_patch >= 0);
  assert(col_patch >= 0);

  int row_idx = (((rowm * amr->im[2] + rowz - amr->ib[2]) *
		   amr->im[1] + rowy - amr->ib[1]) *
		 amr->im[0] + rowx - amr->ib[0]);
  int col_idx = (((colm * amr->im[2] + colz - amr->ib[2]) *
		   amr->im[1] + coly - amr->ib[1]) *
		 amr->im[0] + colx - amr->ib[0]);
  
  if (amr->nr_rows == 0 || amr->rows[amr->nr_rows - 1].idx != row_idx) {
    // start new row
    assert(amr->nr_rows < max_rows);
    amr->rows[amr->nr_rows].patch = row_patch;
    amr->rows[amr->nr_rows].idx = row_idx;
    amr->rows[amr->nr_rows].first_entry = amr->nr_entries;
    amr->nr_rows++;
  }

  assert(amr->nr_entries < max_entries);
  amr->entries[amr->nr_entries].patch = col_patch;
  amr->entries[amr->nr_entries].idx = col_idx;
  amr->entries[amr->nr_entries].val = val;
  amr->nr_entries++;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_assemble

static void
mrc_ddc_amr_assemble(struct mrc_ddc_amr *amr)
{
  amr->rows[amr->nr_rows].first_entry = amr->nr_entries;
  mprintf("nr_rows %d nr_entries %d\n", amr->nr_rows, amr->nr_entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_apply

static void
mrc_ddc_amr_apply(struct mrc_ddc_amr *amr, struct mrc_m3 *fld)
{
  for (int row = 0; row < amr->nr_rows; row++) {
    struct mrc_m3_patch *fldp_row = mrc_m3_patch_get(fld, amr->rows[row].patch);
    int row_idx = amr->rows[row].idx;
    float sum = 0.;
    for (int entry = amr->rows[row].first_entry; entry < amr->rows[row + 1].first_entry; entry++) {
      struct mrc_m3_patch *fldp_col = mrc_m3_patch_get(fld, amr->entries[entry].patch);
      int col_idx = amr->entries[entry].idx;
      float val = amr->entries[entry].val;
      sum += val * fldp_col->arr[col_idx];
    }
    fldp_row->arr[row_idx] = sum;
  }
}

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

struct mrc_ddc_amr_stencil {
  int dx[3];
  float val;
};
  
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
set_ddc_same(struct mrc_ddc_amr *amr, int m, int bnd, int ext[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    int p_nei;
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  mrc_domain_get_neighbor_patch_same(amr->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
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
	    for (int iz = 0; iz < len[2]; iz++) {
	      for (int iy = 0; iy < len[1]; iy++) {
		for (int ix = 0; ix < len[0]; ix++) {
		  mrc_ddc_amr_add_value(amr,
					p,     m, ix + off_to[0]  , iy + off_to[1]  , iz + off_to[2],
					p_nei, m, ix + off_from[0], iy + off_from[1], iz + off_from[2],
					1.f);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

static bool
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

  if (p==9 && i[0] == -1 && ii[1] == 9) {
    printf("dir %d %d dirx %d %d\n", dir[0], dir[1], dirx[0], dirx[1]);
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
mrc_ddc_amr_stencil_coarse(struct mrc_ddc_amr *amr, int ext[3],
			   struct mrc_ddc_amr_stencil *stencil, int stencil_len,
			   int m, int p, int i[3])
{
  int p_nei, j[3];
  mrc_domain_find_valid_point_coarse(amr->domain, ext, p,
				     (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
				     &p_nei, j);
  if (p_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil *s = stencil; s < stencil + stencil_len; s++) {
    int jd[3], p_dnei, j_dnei[3];
    for (int d = 0; d < 3; d++) {
      jd[d] = j[d] + s->dx[d];
    }
    mrc_domain_to_valid_point_same(amr->domain, ext, p_nei, jd, &p_dnei, j_dnei);
    mrc_ddc_amr_add_value(amr, p, m, i[0], i[1], i[2],
			  p_dnei, m, j_dnei[0], j_dnei[1], j_dnei[2], s->val);
  }
  return true;
}

static bool
mrc_ddc_amr_stencil_fine(struct mrc_ddc_amr *amr, int ext[3],
			 struct mrc_ddc_amr_stencil *stencil, int stencil_len,
			 int m, int p, int i[3])
{
  int p_nei, j[3];
  mrc_domain_find_valid_point_fine(amr->domain, p, (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &p_nei, j);
  if (p_nei < 0) {
    return false;
  }

  for (struct mrc_ddc_amr_stencil *s = stencil; s < stencil + stencil_len; s++) {
    int id[3];
    for (int d = 0; d < 3; d++) {
      id[d] = 2*i[d] + s->dx[d];
    }
    mrc_domain_find_valid_point_fine(amr->domain, p, id, &p_nei, j);
    mrc_ddc_amr_add_value(amr, p, m, i[0], i[1], i[2], p_nei, m, j[0], j[1], j[2], s->val);
  }
  return true;
}

// ================================================================================

static void
set_ddc_HZ(struct mrc_ddc_amr *amr, int m, int bnd, int ext[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    for (int iz = 0; iz < ldims[2] + 0; iz++) {
      for (int iy = -bnd; iy < ldims[1] + ext[1] + bnd; iy++) {
	for (int ix = -bnd; ix < ldims[0] + ext[0] + bnd; ix++) {
	  if (ix >= 0 && ix < ldims[0] + ext[0] &&
	      iy >= 0 && iy < ldims[1] + ext[1] &&
	      iz >= 0 && iz < ldims[2] + ext[2]) {
	    assert(!mrc_domain_is_ghost(amr->domain, ext, p, (int[]) { ix, iy, iz }));
	    continue;
	  }
	  if (!mrc_domain_is_ghost(amr->domain, ext, p, (int[]) { ix, iy, iz })) {
	    continue;
	  }

	  int j[3], p_nei;
	  mrc_domain_find_valid_point_same(amr->domain, ext,
					   p, (int[]) { ix, iy, iz }, &p_nei, j);
	  if (p_nei >= 0) {
	    assert(!mrc_domain_is_ghost(amr->domain, ext, p_nei, j));
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], 1.f);
	    continue;
	  }
	  
	  mrc_domain_find_valid_point_coarse(amr->domain, ext,
					     p, (int[]) { div_2(ix), div_2(iy), div_2(iz) }, &p_nei, j);
	  if (p_nei >= 0) {
	    assert(m == HZ);
	    int k = 0; // FIXME, 3D
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2]    , .5f);
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2] + k, .5f);
	    continue;
	  }

	  mrc_domain_find_valid_point_fine(amr->domain, p, (int[]) { 2*ix  , 2*iy  , 2*iz }, &p_nei, j);
	  if (p_nei >= 0) {
	    assert(m == HZ);
	    // FIXME, 3D
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], (2.f/8.f) * 1.f);

	    mrc_domain_find_valid_point_fine(amr->domain, p, (int[]) { 2*ix+1, 2*iy  , 2*iz }, &p_nei, j);
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], (2.f/8.f) * 1.f);

	    mrc_domain_find_valid_point_fine(amr->domain, p, (int[]) { 2*ix  , 2*iy+1, 2*iz }, &p_nei, j);
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], (2.f/8.f) * 1.f);

	    mrc_domain_find_valid_point_fine(amr->domain, p, (int[]) { 2*ix+1, 2*iy+1, 2*iz }, &p_nei, j);
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], (2.f/8.f) * 1.f);
	    continue;
	  }
	  // oops, no way to fill this point?
	  MHERE;
	}
      }
    }
  }
}

static void
set_ddc_EY(struct mrc_ddc_amr *amr, int m, int bnd, int ext[3])
{
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    for (int iz = 0; iz < ldims[2] + 0; iz++) {
      for (int iy = -bnd; iy < ldims[1] + ext[1] + bnd; iy++) {
	for (int ix = -bnd; ix < ldims[0] + ext[0] + bnd; ix++) {
	  if (ix >= ext[0] && ix < ldims[0] &&
	      iy >= ext[1] && iy < ldims[1] &&
	      iz >= ext[2] && iz < ldims[2]) {
	    assert(!mrc_domain_is_ghost(amr->domain, ext, p, (int[]) { ix, iy, iz }));
	    continue;
	  }
	  if (!mrc_domain_is_ghost(amr->domain, ext, p, (int[]) { ix, iy, iz })) {
	    continue;
	  }

	  // at this point, we skipped all interior points, so only ghostpoints are left

	  // try to find an interior point corresponding to the current ghostpoint
	  int j[3], p_nei;
	  mrc_domain_find_valid_point_same(amr->domain, ext,
					   p, (int[]) { ix, iy, iz }, &p_nei, j);
	  if (p_nei >= 0) {
	    assert(!mrc_domain_is_ghost(amr->domain, ext, p_nei, j)); // FIXME, -> add_value()
	    mrc_ddc_amr_add_value(amr, p, m, ix, iy, iz, p_nei, m, j[0], j[1], j[2], 1.f);
	    continue;
	  }

	  // try to interpolate from coarse
	  assert(m == EY);
	  int i = ix & 1;
	  struct mrc_ddc_amr_stencil stencil_coarse[2] = {
	    // FIXME, 3D
	    { .dx = { 0, 0, 0 }, .val = .5f },
	    { .dx = { i, 0, 0 }, .val = .5f },
	  };

	  if (mrc_ddc_amr_stencil_coarse(amr, ext, stencil_coarse, 2, m, p, (int[]) { ix, iy, iz })) {
	    continue;
	  }
	      
	  // try to restrict from fine
	  assert(m == EY);
	  struct mrc_ddc_amr_stencil stencil_fine[6] = {
	    // FIXME, 3D
	    { .dx = { -1,  0,  0 }, .val = (1.f/8.f) * 1.f },
	    { .dx = { -1, +1,  0 }, .val = (1.f/8.f) * 1.f },
	    { .dx = {  0,  0,  0 }, .val = (1.f/8.f) * 2.f },
	    { .dx = {  0, +1,  0 }, .val = (1.f/8.f) * 2.f },
	    { .dx = { +1,  0,  0 }, .val = (1.f/8.f) * 1.f },
	    { .dx = { +1, +1,  0 }, .val = (1.f/8.f) * 1.f },
	  };

	  if (mrc_ddc_amr_stencil_fine(amr, ext, stencil_fine, 6, m, p, (int[]) { ix, iy, iz })) {
	    continue;
	  }
	  // oops, no way to fill this point?
	  MHERE;
	}
      }
    }
  }
}

static void
set_ddc_HY_coarse(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 0, 1, 0 };
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(amr->domain, p, &pi);
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int p_nei;
	  mrc_domain_get_neighbor_patch_coarse(amr->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
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
		len[d] = ldims[d] + ext[d];
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
		  int j = (iy + off_to[1]) & 1;
		  mrc_ddc_amr_add_value(amr,
					p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, HY, (ix + off_from[0])/2, (iy + off_from[1])/2    , (iz + off_from[2])/2,
					.5f);
		  mrc_ddc_amr_add_value(amr, 
					p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, HY, (ix + off_from[0])/2, (iy + off_from[1])/2 + j, (iz + off_from[2])/2, 
					.5f);

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
set_ddc_EZ_coarse(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 1, 1, 0 };
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch_info pi;
    mrc_domain_get_local_patch_info(amr->domain, p, &pi);
    int dir[3];
    for (dir[2] = 0; dir[2] <= 0; dir[2]++) { // FIXME
      for (dir[1] = -1; dir[1] <= 1; dir[1]++) {
	for (dir[0] = -1; dir[0] <= 1; dir[0]++) {
	  if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	    continue;
	  }
	  int p_nei;
	  mrc_domain_get_neighbor_patch_coarse(amr->domain, p, dir, &p_nei);
	  if (p_nei >= 0) {
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
		  int i = (ix + off_to[0]) & 1;
		  int j = (iy + off_to[1]) & 1;
		  mrc_ddc_amr_add_value(amr,
					p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EZ, (ix + off_from[0])/2    , (iy + off_from[1])/2    ,iz,
					.25f);
		  mrc_ddc_amr_add_value(amr,
					p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EZ, (ix + off_from[0])/2 + i, (iy + off_from[1])/2    ,iz,
					.25f);
		  mrc_ddc_amr_add_value(amr,
					p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EZ, (ix + off_from[0])/2    , (iy + off_from[1])/2 + j,iz,
					.25f);
		  mrc_ddc_amr_add_value(amr,
					p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EZ, (ix + off_from[0])/2 + i, (iy + off_from[1])/2 + j,iz,
					.25f);
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
set_ddc_HY_fine(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 0, 1, 0 };
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
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
		mrc_domain_get_neighbor_patch_fine(amr->domain, p, dir, off, &p_nei);
		if (p_nei < 0) {
		  continue;
		}

		for (int iz = 0; iz < 1; iz++) {
		  for (int iy = 0; iy < len[1]; iy++) {
		    for (int ix = 0; ix < len[0]; ix++) {
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix   + off_from[0], 2*iy-1 + off_from[1],iz,
					    (2.f/8.f) * .5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix   + off_from[0], 2*iy   + off_from[1],iz,
					    (2.f/8.f) * 1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix   + off_from[0], 2*iy+1 + off_from[1],iz,
					    (2.f/8.f) * .5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix+1 + off_from[0], 2*iy-1 + off_from[1],iz,
					    (2.f/8.f) * .5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix+1 + off_from[0], 2*iy   + off_from[1],iz,
					    (2.f/8.f) * 1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HY, 2*ix+1 + off_from[0], 2*iy+1 + off_from[1],iz,
					    (2.f/8.f) * .5f);
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
set_ddc_EZ_fine(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 1, 1, 0 };
  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(amr->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
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
		mrc_domain_get_neighbor_patch_fine(amr->domain, p, dir, off, &p_nei);
		if (p_nei < 0) {
		  continue;
		}

		for (int iz = 0; iz < 1; iz++) {
		  for (int iy = 0; iy < len[1]; iy++) {
		    for (int ix = 0; ix < len[0]; ix++) {
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix-1 + off_from[0], 2*iy-1 + off_from[1],iz, (2.f/8.f)*.25f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix-1 + off_from[0], 2*iy   + off_from[1],iz, (2.f/8.f)*.5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix-1 + off_from[0], 2*iy+1 + off_from[1],iz, (2.f/8.f)*.25f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix   + off_from[0], 2*iy-1 + off_from[1],iz, (2.f/8.f)*.5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix   + off_from[0], 2*iy   + off_from[1],iz, (2.f/8.f)*1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix   + off_from[0], 2*iy+1 + off_from[1],iz, (2.f/8.f)*.5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix+1 + off_from[0], 2*iy-1 + off_from[1],iz, (2.f/8.f)*.25f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix+1 + off_from[0], 2*iy   + off_from[1],iz, (2.f/8.f)*.5f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EZ, 2*ix+1 + off_from[0], 2*iy+1 + off_from[1],iz, (2.f/8.f)*.25f);
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

#define F3 MRC_M3

// ======================================================================

static struct mrc_ddc_amr *
make_ddc_H(struct mrc_domain *domain, int sw)
{
  static struct mrc_ddc_amr _ddc;
  struct mrc_ddc_amr *ddc = &_ddc;
  ddc->domain = domain;
  ddc->sw = sw;
  mrc_ddc_amr_setup(ddc);
  set_ddc_same(ddc, HY, 2, (int[]) { 0, 1, 0 });
  set_ddc_HY_coarse(ddc); // not needed for fdtd
  set_ddc_HY_fine(ddc);
  set_ddc_HZ(ddc, HZ, 2, (int[]) { 0, 0, 1 });
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

static struct mrc_ddc_amr *
make_ddc_E(struct mrc_domain *domain, int sw)
{
  static struct mrc_ddc_amr _ddc;
  struct mrc_ddc_amr *ddc = &_ddc;
  ddc->domain = domain;
  ddc->sw = sw;
  mrc_ddc_amr_setup(ddc);
  set_ddc_EY(ddc, EY, 2, (int[]) { 1, 0, 1 });
  set_ddc_same(ddc, EZ, 2, (int[]) { 1, 1, 0 });
  set_ddc_EZ_coarse(ddc);
  set_ddc_EZ_fine(ddc); // not needed for fdtd
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

static void __unused
find_ghosts(struct mrc_domain *domain, struct mrc_m3 *fld, int m,
	    int ext[3], int bnd)
{
  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);
  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    for (int iz = 0; iz < ldims[2] + 0; iz++) {
      for (int iy = -bnd; iy < ldims[1] + ext[1] + bnd; iy++) {
	for (int ix = -bnd; ix < ldims[0] + ext[0] + bnd; ix++) {
	  bool is_ghost = mrc_domain_is_ghost(domain, ext, p, (int[]) { ix, iy, iz });
	  if (!is_ghost) {
	    MRC_M3(fldp, m, ix,iy,iz) = 1.;
	  } else {
	    MRC_M3(fldp, m, ix,iy,iz) = 1./0.;
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// step_fdtd

static void
step_fdtd(struct mrc_m3 *fld, struct mrc_ddc_amr *ddc_E, struct mrc_ddc_amr *ddc_H)
{
  struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
#ifdef AMR
  float dt = 1. / 64;
#else
  float dt = 1. / 16;
#endif

  mrc_ddc_amr_apply(ddc_H, fld);

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, EX, ix,iy,iz) +=
      	cny * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix,iy-1,iz)) -
      	cnz * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix,iy,iz-1));

      F3(fldp, EY, ix,iy,iz) +=
	cnz * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy,iz-1)) -
	cnx * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix-1,iy,iz));

      F3(fldp, EZ, ix,iy,iz) +=
      	cnx * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix-1,iy,iz)) -
      	cny * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy-1,iz));
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  mrc_ddc_amr_apply(ddc_E, fld);

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, HX, ix,iy,iz) -=
	cny * (F3(fldp, EZ, ix,iy+1,iz) - F3(fldp, EZ, ix,iy,iz)) -
	cnz * (F3(fldp, EY, ix,iy,iz+1) - F3(fldp, EY, ix,iy,iz));
      
      F3(fldp, HY, ix,iy,iz) -=
	cnz * (F3(fldp, EX, ix,iy,iz+1) - F3(fldp, EX, ix,iy,iz)) -
	cnx * (F3(fldp, EZ, ix+1,iy,iz) - F3(fldp, EZ, ix,iy,iz));
      
      F3(fldp, HZ, ix,iy,iz) -=
	cnx * (F3(fldp, EY, ix+1,iy,iz) - F3(fldp, EY, ix,iy,iz)) -
	cny * (F3(fldp, EX, ix,iy+1,iz) - F3(fldp, EX, ix,iy,iz));
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, HX, ix,iy,iz) -=
	cny * (F3(fldp, EZ, ix,iy+1,iz) - F3(fldp, EZ, ix,iy,iz)) -
	cnz * (F3(fldp, EY, ix,iy,iz+1) - F3(fldp, EY, ix,iy,iz));
      
      F3(fldp, HY, ix,iy,iz) -=
	cnz * (F3(fldp, EX, ix,iy,iz+1) - F3(fldp, EX, ix,iy,iz)) -
	cnx * (F3(fldp, EZ, ix+1,iy,iz) - F3(fldp, EZ, ix,iy,iz));
      
      F3(fldp, HZ, ix,iy,iz) -=
	cnx * (F3(fldp, EY, ix+1,iy,iz) - F3(fldp, EY, ix,iy,iz)) -
	cny * (F3(fldp, EX, ix,iy+1,iz) - F3(fldp, EX, ix,iy,iz));
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  mrc_ddc_amr_apply(ddc_H, fld);

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    //    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = 0.;//.5 * dt / dy;
    float cnz = 0.;
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      F3(fldp, EX, ix,iy,iz) +=
	cny * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix,iy-1,iz)) -
	cnz * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix,iy,iz-1));
      
      F3(fldp, EY, ix,iy,iz) +=
	cnz * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy,iz-1)) -
	cnx * (F3(fldp, HZ, ix,iy,iz) - F3(fldp, HZ, ix-1,iy,iz));
      
      F3(fldp, EZ, ix,iy,iz) +=
	cnx * (F3(fldp, HY, ix,iy,iz) - F3(fldp, HY, ix-1,iy,iz)) -
	cny * (F3(fldp, HX, ix,iy,iz) - F3(fldp, HX, ix,iy-1,iz));
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

}

static void __unused
set_domain_0(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 1, 1, 0 });
}

static void __unused
set_domain_1(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 1, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

static void __unused
set_domain_2(struct mrc_domain *domain)
{
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
}

static void __unused
set_domain_3(struct mrc_domain *domain)
{
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
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 5, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 5, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

static void __unused
set_domain_4(struct mrc_domain *domain)
{
  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 0, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 0, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 1, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 1, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 4, 3, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 5, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 1, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 2, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 4, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 2, 5, 0 });
  mrc_domain_add_patch(domain, 3, (int [3]) { 3, 5, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 2, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 2, 0 });

  mrc_domain_add_patch(domain, 2, (int [3]) { 0, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 1, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 2, 3, 0 });
  mrc_domain_add_patch(domain, 2, (int [3]) { 3, 3, 0 });
}

float
func1(float x, float y)
{
  float kx = 2. * M_PI;
  return sin(.5 + kx * x);
}

float
func2(float x, float y)
{
  float kx = 2. * M_PI, ky = 2. * M_PI;
  return sin(.5 + kx * x) * cos(.5 + ky * y);
}

float (*func)(float, float) = func2;

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
  set_domain_3(domain);

  mrc_domain_setup(domain);
  mrc_domain_plot(domain);

  // create and fill a field

  struct mrc_m3 *fld = mrc_domain_m3_create(domain);
  mrc_m3_set_name(fld, "fld");
  mrc_m3_set_param_int(fld, "nr_comps", NR_COMPS);
  mrc_m3_set_param_int(fld, "sw", 3);
  mrc_m3_set_from_options(fld);
  mrc_m3_setup(fld);
  mrc_m3_set_comp_name(fld, EX, "EX");
  mrc_m3_set_comp_name(fld, EY, "EY");
  mrc_m3_set_comp_name(fld, EZ, "EZ");
  mrc_m3_set_comp_name(fld, HX, "HX");
  mrc_m3_set_comp_name(fld, HY, "HY");
  mrc_m3_set_comp_name(fld, HZ, "HZ");

  int ldims[3];
  mrc_domain_get_param_int3(fld->domain, "m", ldims);

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);

#if 1
    mrc_m3_foreach(fldp, ix,iy,iz, 3, 3) {
      MRC_M3(fldp, EY, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, HZ, ix,iy,iz) = 1.f / 0.f;
    } mrc_m3_foreach_end;
#endif
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      float x_cc = MRC_MCRDX(crds, ix);
      float y_cc = MRC_MCRDY(crds, iy);
      float x_nc = .5f * (MRC_MCRDX(crds, ix-1) + MRC_MCRDX(crds, ix));
      float y_nc = .5f * (MRC_MCRDY(crds, iy-1) + MRC_MCRDY(crds, iy));
      MRC_M3(fldp, EZ, ix,iy,iz) = func(x_nc, y_nc);
      if (ix < ldims[0]) {
	MRC_M3(fldp, HY, ix,iy,iz) = func(x_cc, y_nc);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 1, 0, 1 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, EY, ix,iy,iz) = func(x_nc, y_cc);
      }

      if (ix < ldims[0] && iy < ldims[1]) {
	MRC_M3(fldp, HZ, ix,iy,iz) = func(x_cc, y_cc);
      }
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  struct mrc_ddc_amr *ddc_E = make_ddc_E(domain, fld->sw);
  struct mrc_ddc_amr *ddc_H = make_ddc_H(domain, fld->sw);

  // write field to disk

  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_type(io, "xdmf2");
  mrc_io_set_param_int(io, "sw", 3);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);

  mrc_io_open(io, "w", 0, 0);
  mrc_m3_write(fld, io);
  mrc_io_close(io);

  mrc_ddc_amr_apply(ddc_E, fld);
  mrc_ddc_amr_apply(ddc_H, fld);
#if 0
  find_ghosts(domain, fld, EY, (int[]) { 1, 0, 1 }, 2);
  find_ghosts(domain, fld, EZ, (int[]) { 1, 1, 0 }, 2);
#endif

  mrc_io_open(io, "w", 1, 1);
  mrc_m3_write(fld, io);
  mrc_io_close(io);

  for (int n = 0; n <= 32; n++) {
    mrc_io_open(io, "w", n+2, n+2);
    mrc_m3_write(fld, io);
    mrc_io_close(io);

    step_fdtd(fld, ddc_E, ddc_H);
  }

  mrc_io_destroy(io);

  mrc_ddc_amr_destroy(ddc_E);
  mrc_ddc_amr_destroy(ddc_H);

  mrc_m3_destroy(fld);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
