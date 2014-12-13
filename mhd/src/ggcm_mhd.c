
#include "ggcm_mhd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_crds_private.h"
#include "ggcm_mhd_crds_gen.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_ic.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_profile.h>

#include <assert.h>
#include <string.h>

#define ggcm_mhd_ops(mhd) ((struct ggcm_mhd_ops *) mhd->obj.ops)

static const char *fldname[_NR_FLDS] = {
  [_RR1 ] = "rr1",
  [_RV1X] = "rv1x",
  [_RV1Y] = "rv1y",
  [_RV1Z] = "rv1z",
  [_UU1 ] = "uu1",
  [_B1X ] = "b1x",
  [_B1Y ] = "b1y",
  [_B1Z ] = "b1z",

  [_RR2 ] = "rr2",
  [_RV2X] = "rv2x",
  [_RV2Y] = "rv2y",
  [_RV2Z] = "rv2z",
  [_UU2 ] = "uu2",
  [_B2X ] = "b2x",
  [_B2Y ] = "b2y",
  [_B2Z ] = "b2z",

  [_YMASK] = "ymask",
  [_ZMASK] = "zmask",
  [_CMSV ] = "cmsv",

  [_RR  ] = "rr",
  [_PP  ] = "pp",
  [_VX  ] = "vx",
  [_VY  ] = "vy",
  [_VZ  ] = "vz",
  [_BX  ] = "bx",
  [_BY  ] = "by",
  [_BZ  ] = "bz",

  [_TMP1] = "tmp1",
  [_TMP2] = "tmp2",
  [_TMP3] = "tmp3",
  [_TMP4] = "tmp4",

  [_FLX ] = "ex",
  [_FLY ] = "ey",
  [_FLZ ] = "ez",

  [_CX  ] = "cx",
  [_CY  ] = "cy",
  [_CZ  ] = "cz",

  [_XTRA1] = "xtra1",
  [_XTRA2] = "xtra2",

  [_RESIS] = "resis",

  [_CURRX] = "currx",
  [_CURRY] = "curry",
  [_CURRZ] = "currz",

  [_RMASK] = "rmask",
};

// ----------------------------------------------------------------------
// ggcm_mhd methods

static void
_ggcm_mhd_create(struct ggcm_mhd *mhd)
{
  mrc_domain_set_type(mhd->domain, "simple");

  ggcm_mhd_crds_set_param_obj(mhd->crds, "domain", mhd->domain);
  ggcm_mhd_step_set_param_obj(mhd->step, "mhd", mhd);
  ggcm_mhd_diag_set_param_obj(mhd->diag, "mhd", mhd);
  ggcm_mhd_bnd_set_param_obj(mhd->bnd, "mhd", mhd);
  ggcm_mhd_ic_set_param_obj(mhd->ic, "mhd", mhd);

  mrc_fld_set_name(mhd->fld, "ggcm_mhd_fld");
  mrc_fld_set_param_obj(mhd->fld, "domain", mhd->domain);
  mrc_fld_set_param_int(mhd->fld, "nr_spatial_dims", 3);
  mrc_fld_dict_add_obj(mhd->fld, "mhd", mhd);
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_state
//
// update C ggcm_mhd state from Fortran common blocks

void
ggcm_mhd_get_state(struct ggcm_mhd *mhd)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->get_state) {
    ops->get_state(mhd);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_set_state
//
// updated Fortran common blocks from C ggcm_mhd state

void
ggcm_mhd_set_state(struct ggcm_mhd *mhd)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->set_state) {
    ops->set_state(mhd);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup_internal

static void
ggcm_mhd_setup_internal(struct ggcm_mhd *mhd)
{
  const int *ghost_dims = mrc_fld_ghost_dims(mhd->fld);
  const int *dims = mrc_fld_dims(mhd->fld);
  for (int d = 0; d < 3; d++) {
    // local domain size
    mhd->im[d] = dims[d];
    // local domain size incl ghost points
    mhd->img[d] = ghost_dims[d];
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_read

static void
_ggcm_mhd_read(struct ggcm_mhd *mhd, struct mrc_io *io)
{
  ggcm_mhd_read_member_objs(mhd, io);

  ggcm_mhd_setup_internal(mhd);
}

static void
ggcm_mhd_setup_amr_domain(struct ggcm_mhd *mhd)
{
  if (mhd->amr == 1) {
    mrc_domain_add_patch(mhd->domain, 0, (int [3]) { 0, 0, 0 });
  } else if (mhd->amr == 2) {
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 1, 0 });
  } else if (mhd->amr == 3) {
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 3, 0 });
  } else if (mhd->amr == 4) {
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 0, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 4, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 5, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 4, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 5, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 1, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 4, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 4, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 5, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 5, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 2, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 3, 0 });
  } else if (mhd->amr == 5) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 2; j++) {
	mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, 0 });
      }
    }
    for (int i = 0; i < 8; i++) {
      for (int j = 4; j < 8; j++) {
	mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, 0 });
      }
    }
  } else if (mhd->amr == 6) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 3; j++) {
	mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, 0 });
      }
    }
    for (int i = 0; i < 8; i++) {
      for (int j = 6; j < 8; j++) {
	mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, 0 });
      }
    }
  } else if (mhd->amr == 7) {
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 0, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 4, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 5, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 4, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 5, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 1, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 4, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 4, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 2, 5, 0 });
    mrc_domain_add_patch(mhd->domain, 3, (int [3]) { 3, 5, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 2, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 2, 0 });
    
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 0, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 1, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 2, 3, 0 });
    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { 3, 3, 0 });
  } else {
    assert(0);
  }
}

static void
_ggcm_mhd_setup(struct ggcm_mhd *mhd)
{
  ggcm_mhd_step_setup_flds(mhd->step);
  for (int m = 0; m < mrc_fld_nr_comps(mhd->fld); m++) {
    mrc_fld_set_comp_name(mhd->fld, m, fldname[m]);
  }

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  // only set the sw on the domain's crds if they're not already set
  if (1||crds->sw == 0) {
    mrc_crds_set_param_int(crds, "sw", mhd->fld->_nr_ghosts);
  }

  if (mhd->amr > 0) {
    ggcm_mhd_setup_amr_domain(mhd);
  }

  ggcm_mhd_setup_member_objs(mhd);
  ggcm_mhd_setup_internal(mhd);

  if (mhd->amr > 0) {
    // FIXME, all leaked
    mhd->ddc_amr_cc = ggcm_mhd_create_amr_ddc(mhd);
    mhd->ddc_amr_E  = ggcm_mhd_create_amr_ddc_E(mhd);
    mhd->ddc_amr_flux_x = ggcm_mhd_create_amr_ddc_flux_x(mhd);
    mhd->ddc_amr_flux_y = ggcm_mhd_create_amr_ddc_flux_y(mhd);
  }
}

// FIXME
#include <mrc_fld_as_double.h>
void mrc_domain_get_neighbor_patch_same(struct mrc_domain *domain, int p,
					int dx[3], int *p_nei);
void mrc_domain_get_neighbor_patch_fine(struct mrc_domain *domain, int gp,
					int dir[3], int off[3], int *gp_nei);
void mrc_domain_get_neighbor_patch_coarse(struct mrc_domain *domain, int gp,
					  int dx[3], int *gp_nei);
void mrc_domain_find_valid_point_coarse(struct mrc_domain *domain, int ext[3],
					int gp, int i[3], int *gp_nei, int j[3]);
void mrc_domain_find_valid_point_same(struct mrc_domain *domain, int ext[3], int gp, int i[3],
				      int *gp_nei, int j[3]);
void mrc_domain_find_valid_point_fine(struct mrc_domain *domain, int ext[3], int gp, int i[3],
				      int *gp_nei, int j[3]);

static inline int
div_2(int i)
{
  // divide by 2, but always round down
  return i >> 1;
  return (i + 10) / 2 - 5;
}

static bool
is_ghost_b(struct mrc_domain *domain, int ext[3], int gp, int i[3])
{
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
  // do we border a coarse domain? (then it's not a ghost point)
  for (dd[2] = 0; dd[2] >= 0; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	int gp_nei;
	mrc_domain_get_neighbor_patch_coarse(domain, gp, dd, &gp_nei);
	if (gp_nei >= 0) {
	  return false;
	}
      }
    }
  }

  // do we border a fine domain? (then it's a ghost point here on the coarse)
  for (dd[2] = 0; dd[2] >= 0; dd[2]--) {
    for (dd[1] = dir[1]; dd[1] >= dir[1] - dirx[1]; dd[1]--) {
      for (dd[0] = dir[0]; dd[0] >= dir[0] - dirx[0]; dd[0]--) {
	if (dd[0] == 0 && dd[1] == 0 && dd[2] == 0) {
	  continue;
	}
	int gp_nei;
	mrc_domain_get_neighbor_patch_fine(domain, gp, dd, (int[3]) { 0,0,0 }, &gp_nei);
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
	/* if (gp == 2) { */
	/*   mprintf("gp %d i %d:%d:%d dd %d:%d:%d gp_nei %d\n", gp, i[0],i[1],i[2], dd[0], dd[1], dd[2], gp_nei); */
	/* } */
	if (gp_nei >= 0) {
	  return gp != gp_nei;
	}
      }
    }
  }
  return true;
}

static mrc_fld_data_t
fill_ghosts_b_one(struct mrc_domain *domain, struct mrc_fld *fld,
		  int m, int i[3], int gp)
{
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int ext[3] = {};
  ext[m] = (gdims[m] > 1);

  // try to find an interior point on the same level corresponding to the current ghostpoint
  int j[3], gp_nei;
  mrc_domain_find_valid_point_same(domain, ext, gp, i, &gp_nei, j);

  if (gp_nei >= 0) { // && gp_nei != gp) {
    /* mprintf("AAA gp %d i %d:%d:%d gp_nei %d j %d:%d:%d (%g)\n", */
    /* 	      gp, i[0], i[1], i[2], gp_nei, j[0], j[1], j[2], */
    /* 	      M3(fld, m, j[0],j[1],j[2], gp_nei)); */
    //mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, 1.f);
  return M3(fld, BX+m, j[0],j[1],j[2], gp_nei);
  }
  
  // is there a corresponding point in a fine domain (also on this face)?
  // then restrict this face to get replace coarse value
  mrc_domain_find_valid_point_fine(domain, ext, gp, (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
  if (gp_nei >= 0) {
    mrc_fld_data_t fact = (m != 0 ? .5f : 1.f) * (m != 1 ? .5f : 1.f);
    mrc_fld_data_t val = 0.f;
    for (int dy = 0; dy <= (m != 1); dy++) {
      for (int dx = 0; dx <= (m != 0); dx++) {
	val += fact * M3(fld, BX+m, j[0]+dx,j[1]+dy,j[2], gp_nei);
      }
    }
    return val;
  }

  // is there a corresponding coarse level underneath?
  mrc_domain_find_valid_point_coarse(domain, ext, gp,
				     (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
				     &gp_nei, j);
  if (gp_nei >= 0) {
    if (i[m] % 2 == 0) {
      return M3(fld, BX+m, j[0],j[1],j[2], gp_nei);
    } else {
      if (m == 0) {
	return .5f * (fill_ghosts_b_one(domain, fld, m, (int [3]) { i[0] - 1, i[1], i[2] }, gp) +
		      fill_ghosts_b_one(domain, fld, m, (int [3]) { i[0] + 1, i[1], i[2] }, gp));
      } else if (m == 1) {
	return .5f * (fill_ghosts_b_one(domain, fld, m, (int [3]) { i[0], i[1] - 1, i[2] }, gp) +
		      fill_ghosts_b_one(domain, fld, m, (int [3]) { i[0], i[1] + 1, i[2] }, gp));
      } else {
	assert(0);
      }
    }
  }

  mprintf("XXX gp %d i %d %d %d\n", gp, i[0], i[1], i[2]);

  // return 2.f;
  assert(0);
}

static void
fill_ghosts_b(struct ggcm_mhd *mhd, struct mrc_fld *fld)
{
  int bnd = 3;
  struct mrc_domain *domain = fld->_domain;

  int ldims[3], gdims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  mrc_domain_get_global_dims(domain, gdims);

  int sw[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = (gdims[d] == 1) ? 0 : bnd;
  }

  for (int m = 0; m < 3; m++) {
    int ext[3] = { 0, 0, 0 };
    ext[m] = (gdims[m] > 1);

    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(domain, p, &info);
      int gp = info.global_patch;
      int i[3];
      for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
	for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	  for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	    // skip points which are definitely in the interior
	    if (i[0] >= ext[0] && i[0] < ldims[0] &&
		i[1] >= ext[1] && i[1] < ldims[1] &&
		i[2] >= ext[2] && i[2] < ldims[2]) {
	      assert(!is_ghost_b(domain, ext, gp, i));
	      continue;
	    }

	    // FIXME, should be unnecessary in the end
	    if (!is_ghost_b(domain, ext, gp, i)) {
	      continue;
	    }

	    // at this point, we skipped all interior points, so only ghostpoints are left
	    mrc_fld_data_t val = fill_ghosts_b_one(domain, fld, m, i, gp);
	    M3(fld, BX+m, i[0],i[1],i[2], gp) = val;
	  }
	}
      }
    }
  }





  return;
#if 0
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	continue;
      }

      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(mhd->domain, p, &info);
      int gp = info.global_patch, *ldims = info.ldims;

      // low side
      int gp_nei;
      int dir[3] = {}, off[3] = {};
      dir[d] = -1;
#if 0
      // if there's a neighbor patch at the same refinement level,
      // the fluxes better be equal already, this for debugging / making sure
      mrc_domain_get_neighbor_patch_same(mhd->domain, gp, dir, &gp_nei);

      if (gp_nei >= 0) {
	/* mprintf("gp %d d %d gp_nei %d\n", gp, d, gp_nei); */
	int p_nei = gp_nei; // FIXME, serial only
	if (d == 1) {
	  for (int iz = 0; iz < ldims[2]; iz++) {
	    for (int ix = 0; ix < ldims[0]; ix++) {
	      mprintf("BY[%d,0,%d] = %g // %g\n", ix, iz,
		      M3(fld, BX + d, ix,0,iz, p),
		      M3(fld, BX + d, ix,ldims[1],iz, p_nei));
	    }
	  }
	}
      }
#endif

#if 0
      off[d] = 1; // for low side
      mrc_domain_get_neighbor_patch_fine(mhd->domain, gp, dir, off, &gp_nei);

      if (gp_nei >= 0) {
	if (d == 1) {
	  //	  mprintf("low gp %d d %d gp_nei %d\n", gp, d, gp_nei);
	  int offx = (gdims[0] > 1), offz = (gdims[2] > 1);
	  for (off[2] = 0; off[2] <= offz; off[2]++) {
	    for (off[0] = 0; off[0] <= offx; off[0]++) {
	      mrc_domain_get_neighbor_patch_fine(mhd->domain, gp, dir, off, &gp_nei);
	      for (int iz = 0; iz < 1; iz++) {
		for (int ix = 0; ix < ldims[0] / 2; ix++) {
		  int iy = 0, iy_nei = ldims[1];
		  int p_nei = gp_nei; // FIXME, serial
		  mrc_fld_data_t val =
		    .5f * (M3(fld, BX + d, ix*2  ,iy_nei,iz*2, p_nei) +
			   M3(fld, BX + d, ix*2+1,iy_nei,iz*2, p_nei));
		  mprintf("BYl gp %d [%d,%d,%d] = %g // %g\n", gp, ix, iy, iz,
			  M3(fld, BX + d, ix + ldims[0]/2 * off[0],iy,iz + ldims[2]/2 * off[2], p),
			  val);
		  /* M3(fld, BX + d, ix + ldims[0]/2 * off[0],iy,iz + ldims[2]/2 * off[2], p) = val; */
		}
	      }
	    }
	  }
	}
      }
#endif

#if 0
      // high side
      dir[d] = 1;
      off[d] = 0; // for high side
      mrc_domain_get_neighbor_patch_fine(mhd->domain, gp, dir, off, &gp_nei);

      if (gp_nei >= 0) {
	if (d == 1) {
	  //	  mprintf("high gp %d d %d gp_nei %d\n", gp, d, gp_nei);
	  int offx = (gdims[0] > 1), offz = (gdims[2] > 1);
	  for (off[2] = 0; off[2] <= offz; off[2]++) {
	    for (off[0] = 0; off[0] <= offx; off[0]++) {
	      mrc_domain_get_neighbor_patch_fine(mhd->domain, gp, dir, off, &gp_nei);
	      for (int iz = 0; iz < 1; iz++) {
		for (int ix = 0; ix < ldims[0] / 2; ix++) {
		  int iy = ldims[1], iy_nei = 0;
		  int p_nei = gp_nei; // FIXME, serial
		  mrc_fld_data_t val =
		    .5f * (M3(fld, BX + d, ix*2  ,iy_nei,iz*2, p_nei) +
			   M3(fld, BX + d, ix*2+1,iy_nei,iz*2, p_nei));
		  mprintf("BYh gp %d [%d,%d,%d] = %g // %g\n", gp, ix, iy, iz,
			  M3(fld, BX + d, ix + ldims[0]/2 * off[0],iy,iz + ldims[2]/2 * off[2], p),
			  val);
		  /* M3(fld, BX + d, ix + ldims[0]/2 * off[0],iy,iz + ldims[2]/2 * off[2], p) = val; */
		}
	      }
	    }
	  }
	}
      }
#endif

      // low side
      dir[d] = -1;
      mrc_domain_get_neighbor_patch_coarse(mhd->domain, gp, dir, &gp_nei);

      if (gp_nei >= 0) {
	for (int iy = -4; iy < 0; iy++) {
	  for (int ix = -4; ix < ldims[0] + 4; ix++) {
	    M3(fld, BX, ix,iy,0, gp) = 1e20;
	    M3(fld, BY, ix,iy,0, gp) = 1e20;
	    M3(fld, BZ, ix,iy,0, gp) = 1e20;
	  }
	}
      }

      // high side
      dir[d] = 1;
      mrc_domain_get_neighbor_patch_coarse(mhd->domain, gp, dir, &gp_nei);

      if (gp_nei >= 0) {
	for (int iy = ldims[1]; iy < ldims[1] + 4; iy++) {
	  for (int ix = -4; ix < ldims[0] + 4; ix++) {
	    M3(fld, BX, ix,iy,0, gp) = 1e20;
	    if (iy > ldims[1]) {
	      M3(fld, BY, ix,iy,0, gp) = 1e20;
	    }
	    M3(fld, BZ, ix,iy,0, gp) = 1e20;
	  }
	}
      }

#if 1
      // low side
      dir[d] = -1;
      mrc_domain_get_neighbor_patch_coarse(mhd->domain, gp, dir, &gp_nei);

      if (gp_nei >= 0) {
	if (d == 1) {
	  int ext[3] = { 0, 1, 0 };
	  //mprintf("low gp %d d %d gp_nei %d\n", gp, d, gp_nei);
	  int iz = 0; {
	    for (int iy = -2; iy < 0; iy += 2) {
	      // FIXME, in the double ghosts, there may be fine values
	      for (int ix = -2; ix < ldims[0] + 2; ix += 2) {
		int i[3] = { ix, iy, iz };
		int gp_nei, j[3];
		mrc_domain_find_valid_point_coarse(mhd->domain, ext, gp,
						   (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
						   &gp_nei, j);
		/* mprintf("gp %d i [%d,%d,%d] gp_nei %d j [%d,%d,%d]\n", */
		/* 	gp, i[0], i[1], i[2], gp_nei, j[0], j[1], j[2]); */

		assert(gp_nei >= 0);

		// FIXME, this overwrites values on the boundary, too
		// left BX
		M3(fld, BX, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BX, j[0]  ,j[1],j[2], gp_nei);
		M3(fld, BX, i[0]  ,i[1]+1,i[2], gp) = M3(fld, BX, j[0]  ,j[1],j[2], gp_nei);

		// right BX
		M3(fld, BX, i[0]+2,i[1]  ,i[2], gp) = M3(fld, BX, j[0]+1,j[1],j[2], gp_nei);
		M3(fld, BX, i[0]+2,i[1]+1,i[2], gp) = M3(fld, BX, j[0]+1,j[1],j[2], gp_nei);

		// bottom BY
		M3(fld, BY, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BY, j[0],j[1]  ,j[2], gp_nei);
		M3(fld, BY, i[0]+1,i[1]  ,i[2], gp) = M3(fld, BY, j[0],j[1]  ,j[2], gp_nei);

		// top BY
		/* M3(fld, BY, i[0]  ,i[1]+2,i[2], gp) = M3(fld, BY, j[0],j[1]+1,j[2], gp_nei); */
		/* M3(fld, BY, i[0]+1,i[1]+2,i[2], gp) = M3(fld, BY, j[0],j[1]+1,j[2], gp_nei); */

		// all BZ
		M3(fld, BZ, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]+1,i[1]  ,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]  ,i[1]+1,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]+1,i[1]+1,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);

		// center BX
		for (int dy = 0; dy <= 1; dy++) {
		  M3(fld, BX, i[0]+1,i[1]+dy,i[2], gp) =
		    .5f * (M3(fld, BX, i[0],i[1]+dy,i[2], gp) + M3(fld, BX, i[0]+2,i[1]+dy,i[2], gp));
		}

		// center BY
		for (int dx = 0; dx <= 1; dx++) {
		  M3(fld, BY, i[0]+dx,i[1]+1,i[2], gp) =
		    .5f * (M3(fld, BY, i[0]+dx,i[1],i[2], gp) + M3(fld, BY, i[0]+dx,i[1]+2,i[2], gp));
		}
	      }
	    }
	  }
	}
      }

      // high side
      dir[d] = 1;
      mrc_domain_get_neighbor_patch_coarse(mhd->domain, gp, dir, &gp_nei);

      if (gp_nei >= 0) {
	if (d == 1) {
	  int ext[3] = { 0, 1, 0 };
	  //mprintf("high gp %d d %d gp_nei %d\n", gp, d, gp_nei);
	  int iz = 0; {
	    for (int iy = ldims[1]; iy < ldims[1] + 2; iy += 2) {
	      // FIXME, in the double ghosts, there may be fine values
	      for (int ix = -2; ix < ldims[0] + 2; ix += 2) {
		int i[3] = { ix, iy, iz };
		int gp_nei, j[3];
		mrc_domain_find_valid_point_coarse(mhd->domain, ext, gp,
						   (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
						   &gp_nei, j);
		/* mprintf("gp %d i [%d,%d,%d] gp_nei %d j [%d,%d,%d]\n", */
		/* 	gp, i[0], i[1], i[2], gp_nei, j[0], j[1], j[2]); */

		assert(gp_nei >= 0);

		// FIXME, this overwrites values on the boundary, too
		// left BX
		M3(fld, BX, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BX, j[0]  ,j[1],j[2], gp_nei);
		M3(fld, BX, i[0]  ,i[1]+1,i[2], gp) = M3(fld, BX, j[0]  ,j[1],j[2], gp_nei);

		// right BX
		M3(fld, BX, i[0]+2,i[1]  ,i[2], gp) = M3(fld, BX, j[0]+1,j[1],j[2], gp_nei);
		M3(fld, BX, i[0]+2,i[1]+1,i[2], gp) = M3(fld, BX, j[0]+1,j[1],j[2], gp_nei);

		// bottom BY
		/* M3(fld, BY, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BY, j[0],j[1]  ,j[2], gp_nei); */
		/* M3(fld, BY, i[0]+1,i[1]  ,i[2], gp) = M3(fld, BY, j[0],j[1]  ,j[2], gp_nei); */

		// top BY
		M3(fld, BY, i[0]  ,i[1]+2,i[2], gp) = M3(fld, BY, j[0],j[1]+1,j[2], gp_nei);
		M3(fld, BY, i[0]+1,i[1]+2,i[2], gp) = M3(fld, BY, j[0],j[1]+1,j[2], gp_nei);

		// all BZ
		M3(fld, BZ, i[0]  ,i[1]  ,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]+1,i[1]  ,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]  ,i[1]+1,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);
		M3(fld, BZ, i[0]+1,i[1]+1,i[2], gp) = M3(fld, BZ, j[0],j[1],j[2], gp_nei);

		// center BX
		for (int dy = 0; dy <= 1; dy++) {
		  M3(fld, BX, i[0]+1,i[1]+dy,i[2], gp) =
		    .5f * (M3(fld, BX, i[0],i[1]+dy,i[2], gp) + M3(fld, BX, i[0]+2,i[1]+dy,i[2], gp));
		}

		// center BY
		for (int dx = 0; dx <= 1; dx++) {
		  M3(fld, BY, i[0]+dx,i[1]+1,i[2], gp) =
		    .5f * (M3(fld, BY, i[0]+dx,i[1],i[2], gp) + M3(fld, BY, i[0]+dx,i[1]+2,i[2], gp));
		}
	      }
	    }
	  }
	}
      }
#endif
    }
  }
#endif
}

void
ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m, float bntim)
{
  if (mhd->amr == 0) {
    mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mhd->domain), m, m + 8, fld);
  } else {
    assert(m == 0);
    mrc_ddc_amr_apply(mhd->ddc_amr_cc, fld);
    fill_ghosts_b(mhd, fld);
  }
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd, fld, m, bntim);
}

int
ggcm_mhd_ntot(struct ggcm_mhd *mhd)
{
  const int *ghost_dims = mrc_fld_ghost_dims(mhd->fld);
  return ghost_dims[0] * ghost_dims[1] * ghost_dims[2];
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_cc

void
ggcm_mhd_get_crds_cc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     float crd_cc[3])
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  crd_cc[0] = MRC_MCRDX(crds, ix, p);
  crd_cc[1] = MRC_MCRDY(crds, iy, p);
  crd_cc[2] = MRC_MCRDZ(crds, iz, p);
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_nc

void
ggcm_mhd_get_crds_nc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     float crd_nc[3])
{
  float crd_cc[3], crd_cc_m[3];
  ggcm_mhd_get_crds_cc(mhd, ix  ,iy  ,iz  , p, crd_cc);
  ggcm_mhd_get_crds_cc(mhd, ix-1,iy-1,iz-1, p, crd_cc_m);

  crd_nc[0] = .5f * (crd_cc_m[0] + crd_cc[0]);
  crd_nc[1] = .5f * (crd_cc_m[1] + crd_cc[1]);
  crd_nc[2] = .5f * (crd_cc_m[2] + crd_cc[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_fc

void
ggcm_mhd_get_crds_fc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     int d, float crd_fc[3])
{
  float crd_cc[3];
  float crd_nc[3];
  ggcm_mhd_get_crds_cc(mhd, ix, iy, iz, p, crd_cc);
  ggcm_mhd_get_crds_nc(mhd, ix, iy, iz, p, crd_nc);

  if (d == 0) {
    // Bx located at i, j+.5, k+.5
    crd_fc[0] = crd_nc[0];
    crd_fc[1] = crd_cc[1];
    crd_fc[2] = crd_cc[2];
  } else if (d == 1) {
    // By located at i+.5, j, k+.5
    crd_fc[0] = crd_cc[0];
    crd_fc[1] = crd_nc[1];
    crd_fc[2] = crd_cc[2];
  } else if (d == 2) {
    // Bz located at i+.5, j+.5, k
    crd_fc[0] = crd_cc[0];
    crd_fc[1] = crd_cc[1];
    crd_fc[2] = crd_nc[2];
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_ec

void
ggcm_mhd_get_crds_ec(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     int d, float crd_ec[3])
{
  float crd_cc[3], crd_nc[3];
  ggcm_mhd_get_crds_cc(mhd, ix,iy,iz, p, crd_cc);
  ggcm_mhd_get_crds_nc(mhd, ix,iy,iz, p, crd_nc);

  if (d == 0) {
    // Ex located at i+.5, j, k
    crd_ec[0] = crd_cc[0];
    crd_ec[1] = crd_nc[1];
    crd_ec[2] = crd_nc[2];
  } else if (d == 1) {
    // Ey located at i, j+.5, k
    crd_ec[0] = crd_nc[0];
    crd_ec[1] = crd_cc[1];
    crd_ec[2] = crd_nc[2];
  } else if (d == 2) {
    // Ez located at i, j, k+.5
    crd_ec[0] = crd_nc[0];
    crd_ec[1] = crd_nc[1];
    crd_ec[2] = crd_cc[2];
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_3d_fld
//
// FIXME, this should cache the fields, rather than creating/destroying
// all the time

struct mrc_fld *
ggcm_mhd_get_3d_fld(struct ggcm_mhd *mhd, int nr_comps)
{
  struct mrc_fld *f = mrc_fld_create(ggcm_mhd_comm(mhd));
  mrc_fld_set_type(f , mrc_fld_type(mhd->fld));
  mrc_fld_set_param_obj(f, "domain", mhd->fld->_domain);
  mrc_fld_set_param_int(f, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(f, "nr_comps", nr_comps);
  mrc_fld_set_param_int(f, "nr_ghosts", mhd->fld->_nr_ghosts);
  mrc_fld_setup(f);

  return f;
}

// ----------------------------------------------------------------------
// ggcm_mhd_put_3d_fld

void
ggcm_mhd_put_3d_fld(struct ggcm_mhd *mhd, struct mrc_fld *f)
{
  mrc_fld_destroy(f);
}

// ----------------------------------------------------------------------
// ggcm_mhd_default_box
//
// This function can be called in a subclass's ::create() function to
// set defaults for non-GGCM, normalized MHD-in-a-box simulations
//
// TODO: This should probably be the default in the first place

void
ggcm_mhd_default_box(struct ggcm_mhd *mhd)
{
  mhd->par.rrnorm = 1.f;
  mhd->par.ppnorm = 1.f;
  mhd->par.vvnorm = 1.f;
  mhd->par.bbnorm = 1.f;
  mhd->par.ccnorm = 1.f;
  mhd->par.eenorm = 1.f;
  mhd->par.resnorm = 1.f;
  mhd->par.tnorm = 1.f;
  mhd->par.diffco = 0.f;

  ggcm_mhd_set_param_float(mhd, "isphere", 0.);
  ggcm_mhd_set_param_float(mhd, "diffsphere", 0.);
  ggcm_mhd_set_param_float(mhd, "speedlimit", 1e9);

  // default to periodic boundary conditions
  ggcm_mhd_bnd_set_type(mhd->bnd, "none");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  // generate MHD solver grid from mrc_crds
  ggcm_mhd_crds_gen_set_type(mhd->crds->crds_gen, "mrc");
}

// ======================================================================
// ggcm_mhd class

static void
ggcm_mhd_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ops_box);
}

static struct mrc_param_select magdiffu_descr[] = {
  { .val = MAGDIFFU_NL1  , .str = "nl1"     },
  { .val = MAGDIFFU_RES1 , .str = "res1"    },
  { .val = MAGDIFFU_CONST, .str = "const"   },
  {},
};

#define VAR(x) (void *)offsetof(struct ggcm_mhd, x)
static struct param ggcm_mhd_descr[] = {
  { "gamma"           , VAR(par.gamm)        , PARAM_FLOAT(1.66667f) },
  { "rrmin"           , VAR(par.rrmin)       , PARAM_FLOAT(.1f)      },
  { "bbnorm"          , VAR(par.bbnorm)      , PARAM_FLOAT(30574.f)  },
  { "vvnorm"          , VAR(par.vvnorm)      , PARAM_FLOAT(6692.98f) },
  { "rrnorm"          , VAR(par.rrnorm)      , PARAM_FLOAT(10000.f)  },
  { "ppnorm"          , VAR(par.ppnorm)      , PARAM_FLOAT(7.43866e8)},
  { "ccnorm"          , VAR(par.ccnorm)      , PARAM_FLOAT(3.81885)  },
  { "eenorm"          , VAR(par.eenorm)      , PARAM_FLOAT(204631.f) },
  { "resnorm"         , VAR(par.resnorm)     , PARAM_FLOAT(53.5848e6)},
  { "tnorm"           , VAR(par.tnorm)       , PARAM_FLOAT(.95189935)},
  { "diffconstant"    , VAR(par.diffco)      , PARAM_FLOAT(.03f)     },
  { "diffthreshold"   , VAR(par.diffth)      , PARAM_FLOAT(.75f)     },
  { "diffsphere"      , VAR(par.diffsphere)  , PARAM_FLOAT(6.f)      },
  { "speedlimit"      , VAR(par.speedlimit)  , PARAM_FLOAT(1500.f)   },
  { "thx"             , VAR(par.thx)         , PARAM_FLOAT(.40f)     },
  { "isphere"         , VAR(par.isphere)     , PARAM_FLOAT(3.0f)     },
  { "timelo"          , VAR(par.timelo)      , PARAM_FLOAT(0.f)      },
  { "d_i"             , VAR(par.d_i)         , PARAM_FLOAT(0.f)      },
  { "dtmin"           , VAR(par.dtmin)       , PARAM_FLOAT(.0002f)   },
  { "dbasetime"       , VAR(par.dbasetime)   , PARAM_DOUBLE(0.)      },
  { "modnewstep"      , VAR(par.modnewstep)  , PARAM_INT(1)          },
  { "magdiffu"        , VAR(par.magdiffu)    , PARAM_SELECT(MAGDIFFU_NL1,
							    magdiffu_descr) },
  { "diff_timelo"     , VAR(par.diff_timelo) , PARAM_FLOAT(0.)       },
  { "diff_swbnd"      , VAR(par.diff_swbnd)  , PARAM_FLOAT(-1e30)    },
  { "diff_obnd"       , VAR(par.diff_obnd)   , PARAM_INT(0)          },

  { "monitor_conservation", VAR(par.monitor_conservation), PARAM_BOOL(false)  },
  { "amr"             , VAR(amr)             , PARAM_INT(0)          },

  { "time"            , VAR(time)            , MRC_VAR_FLOAT         },
  { "dt"              , VAR(dt)              , MRC_VAR_FLOAT         },
  { "istep"           , VAR(istep)           , MRC_VAR_INT           },
  { "timla"           , VAR(timla)           , MRC_VAR_FLOAT         },
  { "dacttime"        , VAR(dacttime)        , MRC_VAR_DOUBLE        },

  { "domain"          , VAR(domain)          , MRC_VAR_OBJ(mrc_domain)        },
  { "fld"             , VAR(fld)             , MRC_VAR_OBJ(mrc_fld)           },
  { "crds"            , VAR(crds)            , MRC_VAR_OBJ(ggcm_mhd_crds)     },
  { "step"            , VAR(step)            , MRC_VAR_OBJ(ggcm_mhd_step)     },
  { "diag"            , VAR(diag)            , MRC_VAR_OBJ(ggcm_mhd_diag)     },
  { "bnd"             , VAR(bnd)             , MRC_VAR_OBJ(ggcm_mhd_bnd)      },
  { "ic"              , VAR(ic)              , MRC_VAR_OBJ(ggcm_mhd_ic)       },

  {},
};
#undef VAR

struct mrc_class_ggcm_mhd mrc_class_ggcm_mhd = {
  .name             = "ggcm_mhd",
  .size             = sizeof(struct ggcm_mhd),
  .param_descr      = ggcm_mhd_descr,
  .init             = ggcm_mhd_init,
  .create           = _ggcm_mhd_create,
  .setup            = _ggcm_mhd_setup,
  .read             = _ggcm_mhd_read,
};

// ----------------------------------------------------------------------
// ts_ggcm_mhd_step_calc_rhs
//
// wrapper to be used in a mrc_ts object

void
ts_ggcm_mhd_step_calc_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_fld)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs;
  struct mrc_fld *fld = (struct mrc_fld *) _fld;
  
  mhd->time = time;
  ggcm_mhd_step_calc_rhs(mhd->step, rhs, fld);
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_step_run
//
// wrapper to be used in a mrc_ts object

void
ts_ggcm_mhd_step_run(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *x = (struct mrc_fld *) _x;
  
  mhd->time = mrc_ts_time(ts);
  mhd->dt = mrc_ts_dt(ts);
  mhd->istep = mrc_ts_step_number(ts);

  ggcm_mhd_step_run(mhd->step, x);

  mrc_ts_set_dt(ts, mhd->dt);
}

// ----------------------------------------------------------------------
// ggcm_mhd_main
//
// Helper function that does most of the work of actually running a
// ggcm_mhd based simulation.
// The subclass of ggcm_mhd, and ggcm_mhd_ic are not explicitly set,
// so the default will be used unless overridden on the command line.
// This typically "just works", since the default will be the class you
// registered before calling this function.

int
ggcm_mhd_main(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);
  ggcm_mhd_register();

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  // set up initial condition

  mpi_printf(MPI_COMM_WORLD, "Setting initial condition...\n");
  ggcm_mhd_ic_run(mhd->ic);
  
  // run time integration

  mpi_printf(MPI_COMM_WORLD, "Starting time integration...\n");
  double time_start = MPI_Wtime();

  struct mrc_ts *ts = mrc_ts_create(mrc_domain_comm(mhd->domain));
  if (ggcm_mhd_step_has_calc_rhs(mhd->step)) {
    mrc_ts_set_type(ts, "rk2");
  } else {
    mrc_ts_set_type(ts, "step");
  }
  mrc_ts_set_context(ts, ggcm_mhd_to_mrc_obj(mhd));

  struct mrc_ts_monitor *mon_output =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output, "ggcm");
  mrc_ts_monitor_set_name(mon_output, "mrc_ts_output");
  mrc_ts_add_monitor(ts, mon_output);

  if (mhd->par.monitor_conservation) {
    struct mrc_ts_monitor *mon_conservation =
      mrc_ts_monitor_create(mrc_ts_comm(ts));
    mrc_ts_monitor_set_type(mon_conservation, "conservation");
    mrc_ts_monitor_set_name(mon_conservation, "mrc_ts_conservation");
    mrc_ts_add_monitor(ts, mon_conservation);
  }

  mrc_ts_set_dt(ts, 1e-6);
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(mhd->fld));
  mrc_ts_set_rhs_function(ts, ts_ggcm_mhd_step_calc_rhs, mhd);
  mrc_ts_set_step_function(ts, ts_ggcm_mhd_step_run, mhd);
  mrc_ts_set_from_options(ts);
  mrc_ts_view(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);

  double time_end = MPI_Wtime();

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int gsize = gdims[0] * gdims[1] * gdims[2];

  double cpu_time = time_end - time_start;
  mpi_printf(MPI_COMM_WORLD,"elapsed time = %g sec.\n", cpu_time);
  mpi_printf(MPI_COMM_WORLD,"\ncell-steps / second = %e\n",
	     (double) gsize * mrc_ts_step_number(ts) / cpu_time);

  mrc_ts_destroy(ts);  

  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

