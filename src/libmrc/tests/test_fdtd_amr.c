
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
const int max_rows = 8000;
const int max_entries = 16000;

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
// mrc_ddc_add_setup

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
    fldp_row->arr[row_idx] = 0.;
    for (int entry = amr->rows[row].first_entry; entry < amr->rows[row + 1].first_entry; entry++) {
      struct mrc_m3_patch *fldp_col = mrc_m3_patch_get(fld, amr->entries[entry].patch);
      int col_idx = amr->entries[entry].idx;
      float val = amr->entries[entry].val;
      fldp_row->arr[row_idx] += val * fldp_col->arr[col_idx];
    }
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
set_ddc_HZ_coarse(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 0, 0, 1 };
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
		  int k = 0;
		  mrc_ddc_amr_add_value(amr,
					p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, HZ, (ix + off_from[0])/2, (iy + off_from[1])/2, (iz + off_from[2])/2    ,
					.5f);
		  mrc_ddc_amr_add_value(amr,
					p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, HZ, (ix + off_from[0])/2, (iy + off_from[1])/2, (iz + off_from[2])/2 + k,
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
set_ddc_EY_coarse(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 1, 0, 1 };
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
		  int k = (ix + off_to[0]) & 1;
		  mrc_ddc_amr_add_value(amr,
					p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EY, (ix + off_from[0])/2    , (iy + off_from[1])/2,iz,
					.5f);
		  mrc_ddc_amr_add_value(amr,
					p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					p_nei, EY, (ix + off_from[0])/2 + k, (iy + off_from[1])/2,iz,
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
set_ddc_HZ_fine(struct mrc_ddc_amr *amr)
{
  int bnd = 2;
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
		mrc_domain_get_neighbor_patch_fine(amr->domain, p, dir, off, &p_nei);
		if (p_nei < 0) {
		  continue;
		}

		for (int iz = 0; iz < 1; iz++) {
		  for (int iy = 0; iy < len[1]; iy++) {
		    for (int ix = 0; ix < len[0]; ix++) {
		      mrc_ddc_amr_add_value(amr,
					    p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HZ, 2*ix   + off_from[0], 2*iy   + off_from[1],iz,
					    (2.f/8.f) * 1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HZ, 2*ix+1 + off_from[0], 2*iy   + off_from[1],iz,
					    (2.f/8.f) * 1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HZ, 2*ix   + off_from[0], 2*iy+1 + off_from[1],iz,
					    (2.f/8.f) * 1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     HZ, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, HZ, 2*ix+1 + off_from[0], 2*iy+1 + off_from[1],iz,
					    (2.f/8.f) * 1.f);
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
set_ddc_EY_fine(struct mrc_ddc_amr *amr)
{
  int bnd = 2, ext[3] = { 1, 0, 1 };
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
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix-1 + off_from[0], 2*iy   + off_from[1],iz, (1.f/8.f)*1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix-1 + off_from[0], 2*iy+1 + off_from[1],iz, (1.f/8.f)*1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix   + off_from[0], 2*iy   + off_from[1],iz, (1.f/8.f)*2.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix   + off_from[0], 2*iy+1 + off_from[1],iz, (1.f/8.f)*2.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix-1 + off_from[0], 2*iy   + off_from[1],iz, (1.f/8.f)*1.f);
		      mrc_ddc_amr_add_value(amr,
					    p,     EY, ix + off_to[0], iy + off_to[1], iz + off_to[2],
					    p_nei, EY, 2*ix-1 + off_from[0], 2*iy+1 + off_from[1],iz, (1.f/8.f)*1.f);
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
  set_ddc_same(ddc, HZ, 2, (int[]) { 0, 0, 1 });
  set_ddc_HY_coarse(ddc); // not needed for fdtd
  set_ddc_HZ_coarse(ddc); // not needed for fdtd
  set_ddc_HY_fine(ddc);
  set_ddc_HZ_fine(ddc);
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
  set_ddc_same(ddc, EY, 2, (int[]) { 1, 0, 1 });
  set_ddc_same(ddc, EZ, 2, (int[]) { 1, 1, 0 });
  set_ddc_EY_coarse(ddc);
  set_ddc_EZ_coarse(ddc);
  set_ddc_EY_fine(ddc); // not needed for fdtd
  set_ddc_EZ_fine(ddc); // not needed for fdtd
  mrc_ddc_amr_assemble(ddc);

  return ddc;
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
#if 1
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

  float kx = 2. * M_PI;//, ky = 2. * M_PI;

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
      //      float y_cc = MRC_MCRDY(crds, iy);
      float x_nc = .5f * (MRC_MCRDX(crds, ix-1) + MRC_MCRDX(crds, ix));
      //      float y_nc = .5f * (MRC_MCRDY(crds, iy-1) + MRC_MCRDY(crds, iy));
      MRC_M3(fldp, EZ, ix,iy,iz) = sin(.5+kx * x_nc);// * cos(.5+ky * y_nc);
      if (ix < ldims[0]) {
	MRC_M3(fldp, HY, ix,iy,iz) = sin(.5+kx * x_cc);// * cos(.5+ky * y_nc);
      }
      if (iy < ldims[1]) {
	MRC_M3(fldp, EY, ix,iy,iz) = sin(.5+kx * x_nc);// * cos(.5+ky * y_cc);
	if (ix < ldims[0]) {
	  MRC_M3(fldp, HZ, ix,iy,iz) = sin(.5+kx * x_cc);// * cos(.5+ky * y_cc);
	}
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

  for (int n = 0; n <= 32; n++) {
    mrc_ddc_amr_apply(ddc_E, fld); // for debugging only!
    mrc_ddc_amr_apply(ddc_H, fld); // for debugging only!

    mrc_io_open(io, "w", n, n);
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
