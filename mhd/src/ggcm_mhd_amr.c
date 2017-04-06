
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_step.h" // FIXME

#include <mrc_ddc.h>
#include <mrc_domain.h>

#include <string.h>
#include <math.h>

// FIXME -> header
#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

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

void mrc_ddc_amr_add_diagonal_one(struct mrc_ddc *ddc, int gp, int m, int i[3]);

// ======================================================================
//
// Overall, we maintain three ddc's for AMR operation:
//
// ddc_amr_cc: (FIXME, name)
//   This is for fill_ghosts(), and is responsible for communicating
//   boundary points on the same level (e.g., filling ghost points from
//   a neighboring patch's interior point, but also filling a value on the
//   boundary from a neighboring patch which is the "real" owner of that
//   point.) For AMR, it also takes care of restricting/prolongating the
//   solution to fill ghost points / boundary points.
//
// ddc_amr_flux_*:
//   This is for flux correction right on the boundary between coarse and
//   fine: coarse fluxes are corrected to the aggregated fine values on
//   the same face.
//
// ddc_amr_E:
//   This is for correcting the E field on edges that are on a boundary 
//   between coarse and fine: coarse E fields are to the aggregated fine
//   values on the same edge



// ======================================================================
// correct fluxes
// ======================================================================

static struct mrc_ddc_amr_stencil_entry stencil_fine_flux_x[] = {
  { .dx = {  0,  0,  0 }, .val = .25f },
  { .dx = {  0, +1,  0 }, .val = .25f },
  { .dx = {  0,  0, +1 }, .val = .25f },
  { .dx = {  0, +1, +1 }, .val = .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_flux_y[] = {
  { .dx = {  0,  0,  0 }, .val = .25f },
  { .dx = { +1,  0,  0 }, .val = .25f },
  { .dx = {  0,  0, +1 }, .val = .25f },
  { .dx = { +1,  0, +1 }, .val = .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_flux_z[] = {
  { .dx = {  0,  0,  0 }, .val = .25f },
  { .dx = { +1,  0,  0 }, .val = .25f },
  { .dx = {  0, +1,  0 }, .val = .25f },
  { .dx = { +1, +1,  0 }, .val = .25f },
};

static struct mrc_ddc_amr_stencil stencils_fine_flux[3] = {
  { stencil_fine_flux_x, ARRAY_SIZE(stencil_fine_flux_x) },
  { stencil_fine_flux_y, ARRAY_SIZE(stencil_fine_flux_y) },
  { stencil_fine_flux_z, ARRAY_SIZE(stencil_fine_flux_z) },
};

// ----------------------------------------------------------------------
// ggcm_mhd_create_amr_ddc_flux

struct mrc_ddc *
ggcm_mhd_create_amr_ddc_flux(struct ggcm_mhd *mhd, int d)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_nd->size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  // FIXME!!!
  int n_comp;
  if (strcmp(ggcm_mhd_step_type(mhd->step), "vl") == 0 ||
      strcmp(ggcm_mhd_step_type(mhd->step), "c3_double") == 0) {
    n_comp = 5;
  } else if (strcmp(ggcm_mhd_step_type(mhd->step), "vlct") == 0) {
    n_comp = 8;
  } else {
    assert(0);
  }
  mrc_ddc_set_param_int(ddc, "n_comp", n_comp);
  mrc_ddc_setup(ddc);
  int ext[3] = {};
  ext[d] = 1;
  for (int m = 0; m < n_comp; m++) {
    mrc_ddc_amr_set_by_stencil(ddc, m, (m >= 5) ? 1 : 0, ext, NULL, &stencils_fine_flux[d]);
  }
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// ======================================================================
// fill ghosts
// ======================================================================

static struct mrc_ddc_amr_stencil_entry stencil_coarse_cc[2] = {
  // FIXME, needs some interpolation
  { .dx = { 0, 0, 0 }, .val = 1.f },
};

static struct mrc_ddc_amr_stencil stencils_coarse_cc = {
  stencil_coarse_cc, ARRAY_SIZE(stencil_coarse_cc)
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_cc[] = {
  { .dx = {  0,  0,  0 }, .val = .125f },
  { .dx = { +1,  0,  0 }, .val = .125f },
  { .dx = {  0, +1,  0 }, .val = .125f },
  { .dx = { +1, +1,  0 }, .val = .125f },
  { .dx = {  0,  0, +1 }, .val = .125f },
  { .dx = { +1,  0, +1 }, .val = .125f },
  { .dx = {  0, +1, +1 }, .val = .125f },
  { .dx = { +1, +1, +1 }, .val = .125f },
};

static struct mrc_ddc_amr_stencil stencils_fine_cc = {
  stencil_fine_cc, ARRAY_SIZE(stencil_fine_cc)
};

static inline int
div_2(int i)
{
  // divide by 2, but always round down
  return i >> 1;
}

// ----------------------------------------------------------------------
// fill_ghosts_b_add_one

static void
fill_ghosts_b_add_one(struct mrc_ddc *ddc, struct mrc_domain *domain, int m,
		      int gp0, int i0[3], double scale, int gp, int i[3])
{
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int ext[3] = {};
  ext[m] = (gdims[m] > 1);

  // try to find an interior point on the same level corresponding to the current ghostpoint
  // this also covers the case where the point on the boundary face is a valid value itself
  // (because we border a coarser patch)
  int j[3], gp_nei;
  mrc_domain_find_valid_point_same(domain, ext, gp, i, &gp_nei, j);

  if (gp_nei >= 0) {
    mrc_ddc_amr_add_value(ddc, gp0, BX+m, i0, gp_nei, BX+m, j, scale);
    return;
  }
  
  // is there a corresponding point in a fine domain (also on this face)?
  // then restrict this face to get replacement coarse value
  mrc_domain_find_valid_point_fine(domain, ext, gp, (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
  if (gp_nei >= 0) {
    int nx = (gdims[0] > 1) ? (m != 0) : 0;
    int ny = (gdims[1] > 1) ? (m != 1) : 0;
    int nz = (gdims[2] > 1) ? (m != 2) : 0;
    double fact = (nx ? .5 : 1.) * (ny ? .5 : 1.) * (nz ? .5 : 1.);
    for (int dz = 0; dz <= nz; dz++) {
      for (int dy = 0; dy <= ny; dy++) {
	for (int dx = 0; dx <= nx; dx++) {
	  mrc_ddc_amr_add_value(ddc, gp0, BX+m, i0, gp_nei, BX+m,
				(int[3]) { j[0]+dx,j[1]+dy,j[2]+dz }, fact * scale);
	}
      }
    }
    return;
  }

  // is there a corresponding coarse level underneath?
  mrc_domain_find_valid_point_coarse(domain, ext, gp,
				     (int[]) { div_2(i[0]), div_2(i[1]), div_2(i[2]) },
				     &gp_nei, j);
  if (gp_nei >= 0) {
    if (i[m] % 2 == 0) { // there is a corresponding coarse face that we're on
      // 0th order interpolation
      mrc_ddc_amr_add_value(ddc, gp0, BX+m, i0, gp_nei, BX+m, j, scale);
    } else {
      // we're in the middle between the two underlying coarse faces, so
      // we're basing our value on the already interpolated-to-fine values on those two
      // faces (and in fact, they may be not interpolated values, but actual fine values
      // from a neighboring fine patch)
      if (m == 0) {
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0] - 1, i[1], i[2] });
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0] + 1, i[1], i[2] });
      } else if (m == 1) {
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0], i[1] - 1, i[2] });
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0], i[1] + 1, i[2] });
      } else if (m == 2) {
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0], i[1], i[2] - 1 });
	fill_ghosts_b_add_one(ddc, domain, m, gp0, i, .5 * scale, gp, (int [3]) { i[0], i[1], i[2] + 1 });
      } else {
	assert(0);
      }
    }
    return;
  }

  //mprintf("XXX gp %d i %d %d %d\n", gp, i[0], i[1], i[2]);
  //  assert(0);
}

// ----------------------------------------------------------------------
// ggcm_mhd_create_amr_ddc
//
// for fill_ghosts() when using AMR

struct mrc_ddc *
ggcm_mhd_create_amr_ddc(struct ggcm_mhd *mhd)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_nd->size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  mrc_ddc_set_param_int(ddc, "n_comp", mhd->fld->_nr_comps);
  mrc_ddc_setup(ddc);
  int bnd = mrc_fld_spatial_sw(mhd->fld)[0];

  // cell-center hydro variables
  for (int m = 0; m < 5; m++) {
    mrc_ddc_amr_set_by_stencil(ddc, m, bnd, (int[]) { 0, 0, 0 },
			       &stencils_coarse_cc, &stencils_fine_cc);
  }

  // B field
  int nr_patches, gdims[3];
  mrc_domain_get_patches(mhd->domain, &nr_patches);
  mrc_domain_get_global_dims(mhd->domain, gdims);

  for (int m = 0; m < 3; m++) {
    int ext[3] = { 0, 0, 0 };
    ext[m] = (gdims[m] > 1);

    // when we have bnd = 4 cell-centered ghost cells, we
    // can fill, eg., Bx up to 3 points beyond boundary, but By, Bz up to 4 points
    int sw[3];
    for (int d = 0; d < 3; d++) {
      sw[d] = (ext[d] == 1) ? bnd - 1 : bnd;
      if (gdims[d] == 1) {
	sw[d] = 0;
      }
    }


    for (int p = 0; p < nr_patches; p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(mhd->domain, p, &info);
      int gp = info.global_patch, *ldims = info.ldims;
      
      int i[3];
      for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
	for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	  for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	    // skip points which are definitely in the interior (x)
	    // (not on a boundary face)
	    if (i[0] >= ext[0] && i[0] < ldims[0] &&
		i[1] >= ext[1] && i[1] < ldims[1] &&
		i[2] >= ext[2] && i[2] < ldims[2]) {
	      // E.g., for Bx
	      // +---+---+---+ 
	      // X   x   x   X
	      // +---+---+---+
	      // X   x   x   X
	      // +---+---+---+
	      // X   x   x   X
	      // +---+---+---+ 
	      // --> x
	      mrc_ddc_amr_add_diagonal_one(ddc, gp, BX+m, i);
	      continue;
	    }

	    // at this point, we skipped all interior points, so only values on
	    // and beyond the boundary (ie., ghosts) are left
	    fill_ghosts_b_add_one(ddc, mhd->domain, m, gp, i, 1., gp, i);
	  }
	}
      }
    }
  }
	    
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// ======================================================================
// correct E field
// ======================================================================

// ----------------------------------------------------------------------
// ggcm_mhd_create_amr_ddc_E

struct mrc_ddc *
ggcm_mhd_create_amr_ddc_E(struct ggcm_mhd *mhd)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_nd->size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  mrc_ddc_set_param_int(ddc, "n_comp", 3);
  mrc_ddc_setup(ddc);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int sw[3] = {};

  int nr_patches;
  mrc_domain_get_patches(mhd->domain, &nr_patches);

  for (int m = 0; m < 3; m++) {
    int ext[3];
    for (int d = 0; d < 3; d++) {
      ext[d] = (gdims[d] > 1) ? (m != d) : 0;
    }
    
    for (int p = 0; p < nr_patches; p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(mhd->domain, p, &info);
      int gp = info.global_patch, *ldims = info.ldims;
      
      int i[3];
      for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
	for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	  for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	    if (i[0] >= ext[0] && i[0] < ldims[0] &&
		i[1] >= ext[1] && i[1] < ldims[1] &&
		i[2] >= ext[2] && i[2] < ldims[2]) {
	      // truly interior point "x", ie., not on boundary
	      // X---X---X---X 
	      // |   |   |   |
	      // X---x---x---X
	      // |   |   |   |
	      // X---x---x---X
	      // |   |   |   |
	      // X---X---X---X 
	      mrc_ddc_amr_add_diagonal_one(ddc, gp, m, i);
	      continue;
	    }

	    // now we're only looking at E field values that are on edges on the boundary of
	    // this patch
	    // FIXME, with "c3" double this is really needed though it shouldn't be!!!
	    // to maintain div B = 0
	    
	    int j[3], gp_nei;
	    mrc_domain_find_valid_point_same(mhd->domain, ext, gp, i, &gp_nei, j);
	    if (gp_nei >= 0) {
	      mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, 1.f);
	      continue;
	    }
	    
	    // If we're bordering a fine patch so that this E value is also on the fine grid,
	    // replace this coarse grid value by the one from the fine grid in the same place
	    
	    mrc_domain_find_valid_point_fine(mhd->domain, ext, gp,
					     (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
	    if (gp_nei >= 0) {
	      mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, .5f);
	      mrc_domain_find_valid_point_fine(mhd->domain, ext, gp,
					       (int[]) { 2*i[0] + (m == 0),
							 2*i[1] + (m == 1),
							 2*i[2] + (m == 2) },
					       &gp_nei, j);
	      assert(gp_nei >= 0);
	      mrc_ddc_amr_add_value(ddc, gp, m, i, gp_nei, m, j, .5f);
	      continue;
	    }
	    
	    // FIXME not handled -- should be fixed, though (see FIXME above) not strictly necessary
	    MHERE;
	  }
	}
      }
    }
  }

  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// FIXME, hardcoded field type but for reference / debugging implementation only
#include <mrc_fld_as_double.h>

// ======================================================================
// Reference implementation of B ghost point filling -- now superseded
// by ddc_amr
// ======================================================================

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
    mrc_fld_data_t fact = (m != 0 ? .5f : 1.f) * (m != 1 ? .5f : 1.f) * (m != 2 ? .5f : 1.f);
    mrc_fld_data_t val = 0.f;
    for (int dz = 0; dz <= (m != 2); dz++) {
      for (int dy = 0; dy <= (m != 1); dy++) {
	for (int dx = 0; dx <= (m != 0); dx++) {
	  val += fact * M3(fld, BX+m, j[0]+dx,j[1]+dy,j[2]+dz, gp_nei);
	}
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

  mprintf("XXX m%d gp %d i %d %d %d\n", m, gp, i[0], i[1], i[2]);

  // return 2.f;
  assert(0);
}

// ----------------------------------------------------------------------
// is_ghost_b
//
// determines whether a given B value (on a face) is a ghost value, ie.,
// to be calculated from other (non-ghost) values
// FIXME, does this get double ghost points right? (it may)

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
      dir[1] == 0 && dirx[1] == 0 &&
      dir[2] == 0 && dirx[2] == 0) {
    return false;
  }

  // on the boundary...
  int dd[3];
  // do we border a coarse domain? (then it's not a ghost point)
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
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
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
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

  // is another same level patch in line before us, then it's theirs, and we have
  // a ghost point
  for (dd[2] = dir[2]; dd[2] >= dir[2] - dirx[2]; dd[2]--) {
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

void
ggcm_mhd_amr_fill_ghosts_b(struct ggcm_mhd *mhd, struct mrc_fld *fld)
{
  int bnd = fld->_nr_ghosts - 1;
  struct mrc_domain *domain = mhd->domain;

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
	    /* if (val != M3(fld, BX+m, i[0],i[1],i[2], gp)) { */
	    /*   mprintf("m %d gp %d i %d:%d:%d %g // %g\n", m, gp, i[0],i[1],i[2], val, */
	    /* 	      M3(fld, BX+m, i[0],i[1],i[2], gp)); */
	    /* } */
	    M3(fld, BX+m, i[0],i[1],i[2], gp) = val;
	  }
	}
      }
    }
  }
}


void
correct_E(struct ggcm_mhd *mhd, struct mrc_fld *E)
{
  int bnd = 0;
  struct mrc_domain *domain = mhd->domain;

  int ldims[3], gdims[3];
  mrc_domain_get_param_int3(domain, "m", ldims);
  mrc_domain_get_global_dims(domain, gdims);

  int sw[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = (gdims[d] == 1) ? 0 : bnd;
  }

  for (int m = 0; m < 3; m++) {
    int ext[3];
    for (int d = 0; d < 3; d++) {
      ext[d] = m != d;
    }
    
    for (int p = 0; p < mrc_fld_nr_patches(E); p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(mhd->domain, p, &info);
      int gp = info.global_patch, *ldims = info.ldims;
      
      int i[3];
      for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
	for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	  for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	    if (i[0] >= ext[0] && i[0] < ldims[0] &&
		i[1] >= ext[1] && i[1] < ldims[1] &&
		i[2] >= ext[2] && i[2] < ldims[2]) {
	      // truly interior point "x", ie., not on boundary
	      // X---X---X---X 
	      // |   |   |   |
	      // X---x---x---X
	      // |   |   |   |
	      // X---x---x---X
	      // |   |   |   |
	      // X---X---X---X 
	      continue;
	    }

	    // now we're only looking at E field values that are on edges on the boundary of
	    // this patch
	    // FIXME, with "c3" double this is really needed though it shouldn't be!!!
	    // to maintain div B = 0
	    
	    int j[3], gp_nei;
	    mrc_domain_find_valid_point_same(mhd->domain, ext, gp, i, &gp_nei, j);
	    if (gp_nei >= 0) {
	      assert(M3(E, m, i[0],i[1],i[2], gp) == M3(E, m, j[0],j[1],j[2], gp_nei));
	      continue;
	    }
	    
	    // If we're bordering a fine patch so that this E value is also on the fine grid,
	    // replace this coarse grid value by the one from the fine grid in the same place
	    
	    mrc_domain_find_valid_point_fine(mhd->domain, ext, gp,
					     (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
	    if (gp_nei >= 0) {
	       mrc_fld_data_t val1 = .5f * M3(E, m, j[0],j[1],j[2], gp_nei);
	      int j2[3], gp_nei2;
	      mrc_domain_find_valid_point_fine(mhd->domain, ext, gp,
					       (int[]) { 2*i[0] + (m == 0),
							 2*i[1] + (m == 1),
							 2*i[2] + (m == 2) },
					       &gp_nei2, j2);
	      assert(gp_nei >= 0);
	      mrc_fld_data_t val2 = .5f * M3(E, m, j2[0],j2[1],j2[2], gp_nei2);
	      if (fabs(M3(E, m, i[0],i[1],i[2], gp) - (val1 + val2)) >= 1e-9) {
		mprintf("gp_nei %d j %d:%d:%d is_ghost %d\n", gp_nei, j[0], j[1], j[2],
			mrc_domain_is_ghost(domain, ext, gp_nei, j));
		mprintf("gp %d i %d:%d:%d E %g val %g %g %g j %d %d:%d:%d j2 %d %d:%d:%d\n",
			gp, i[0], i[1], i[2],
			M3(E, m, i[0],i[1],i[2], gp), 2*(val1 + val2), 2*val1, 2*val2,
			gp_nei, j[0], j[1], j[2], gp_nei2, j2[0], j2[1], j2[2]);
		assert(0);
	      }
	      continue;
	    }
	    // FIXME not handled -- should be fixed, though (see FIXME above) not strictly necessary
	    MHERE;
	  }
	}
      }
    }
  }
}

static void
skip_comments(FILE *f) {
  fpos_t preread_pos;
  const int bufflen = 128;
  char line[bufflen];
  int armed = 999;
  int arm_next = 1;

  while (!feof(f)) {
    armed = arm_next;
    fgetpos(f, &preread_pos);
    fgets(line, bufflen, f);

    if (strlen(line) == 0 || line[0] == '\n') {
      continue;
    }
    if (line[strlen(line) - 1] == '\n') {
      arm_next = 1;
    } else {
      arm_next = 0;
    }

    if (armed) {
      if (line[0] != '#') {
	break;
      }
    } else {
      continue;
    }
  }

  if (armed != 999) {
    fsetpos(f, &preread_pos);
  }
}

void
ggcm_mhd_setup_amr_domain(struct ggcm_mhd *mhd)
{
  MPI_Comm mhd_comm = mhd->obj.comm;
  int mhd_rank;  
  MPI_Comm_rank(mhd_comm, &mhd_rank);

  if (mhd->amr == 999) {
    // // FIXME: there should be some feedback if the grid resolution / global
    // // domain size are not the same as the grid file expects
    // struct mrc_patch_info l_info;
    // mrc_domain_get_local_patch_info(mhd->domain, 0, &l_info);

    FILE *f = NULL;
    int ret, mx, my, mz, npatches;
    int level, ix, iy, iz;
    int buf[4];
    double tmp[3];
    int i;
    
    // open the grid file and read the header info
    if (mhd_rank == 0) {
      f = fopen(mhd->amr_grid_file, "r");
      
      if (f == NULL) {
	mprintf("ERROR: tried to make amr grid from '%s', but file DNE.\n",
		mhd->amr_grid_file);
	MPI_Abort(mhd_comm, 1);
      }
    
      skip_comments(f);
      
      ret = fscanf(f, "Level 0 xl: %lg, %lg, %lg\n", &tmp[0], &tmp[1], &tmp[2]);
      if (ret != 3) {MPI_Abort(mhd_comm, 2);}  // malformed file      
      // if (tmp[0] != g_info.xl[0] ||
      //     tmp[1] != g_info.xl[1] || 
      //     tmp[2] != g_info.xl[2]) {
      //   mprintf("WARNING: AMR grid in file expects level 0 to have xl = "
      //           "(%lg, %lg, %lg), but it is actually (%lg, %lg, %lg)\n",
      //           tmp[0], tmp[1], tmp[2],
      //           g_info.xl[0], g_info.xl[1], g_info.xl[2]);
      // }
      ret = fscanf(f, "Level 0 xh: %lg, %lg, %lg\n", &tmp[0], &tmp[1], &tmp[2]);
      if (ret != 3) {MPI_Abort(mhd_comm, 2);}  // malformed file      
      // if (tmp[0] != g_info.xh[0] ||
      //     tmp[1] != g_info.xh[1] ||
      //     tmp[2] != g_info.xh[2]) {
      //   mprintf("WARNING: AMR grid in file expects level 0 to have xh = "
      //           "(%lg, %lg, %lg), but it is actually (%lg, %lg, %lg)\n",
      //           tmp[0], tmp[1], tmp[2], g_info.xh[0], g_info.xh[1], g_info.xh[2]);
      // }    
    
      ret = fscanf(f, "resolution: %d, %d, %d\n", &mx, &my, &mz);
      // mprintf(">> resolution: %d %d %d\n", mx, my, mz);
      if (ret != 3) {MPI_Abort(mhd_comm, 2);}  // malformed file
      // if (mx != l_info.ldims[0] ||
      //     my != l_info.ldims[1] ||
      //     mz != l_info.ldims[2]) {
      //   mprintf("WARNING: AMR grid in file expects patches to have resolution "
      //           "(%d, %d, %d), but it is actually (%d, %d, %d)\n",
      //           mx, my, mz,
      //           l_info.ldims[0], l_info.ldims[1], l_info.ldims[2]);
      // }
      
      ret = fscanf(f, "npatches: %d\n", &npatches);
      // mprintf(">> npatches: %d\n", npatches);
      if (ret != 1) {MPI_Abort(mhd_comm, 3);}  // malformed file
      
      skip_comments(f);
    }
    MPI_Bcast(&npatches, 1, MPI_INT, 0, mhd_comm);
    
    // now add the patches in the file
    for (i = 0; i < npatches + 1; i++) {
      if (mhd_rank == 0) {
	// buf is an array of [level, ix, iy, iz]
	ret = fscanf(f, "%d; %d, %d, %d\n", &buf[0], &buf[1], &buf[2], &buf[3]);
	// mprintf(">> patch: %d; %d %d %d\n", buf[0], buf[1], buf[2], buf[3]);
	if (ret != 4) {
	  buf[0] = -1;  // either the file is malformed, or we're at the end
	}
      }
      MPI_Bcast(buf, 4, MPI_INT, 0, mhd_comm);
      level = buf[0];
      ix = buf[1];
      iy = buf[2];
      iz = buf[3];
      if (level == -1) {
	break;
      }
      mrc_domain_add_patch(mhd->domain, level, (int [3]) {ix, iy, iz}); 
    }
    
    if (mhd_rank == 0) {
      fclose(f);
      if (i != npatches){
	printf("ERROR: amr grid file lied about number of patches\n");
	MPI_Abort(mhd_comm, 4);
      }
    }
  } else if (mhd->amr == 1) {
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
  } else if (mhd->amr == 8) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 4; j++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 9) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 2; j++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
    for (int k = 0; k < 8; k++) {
      for (int i = 0; i < 8; i++) {
	for (int j = 4; j < 8; j++) {
	  mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 10) {
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 0, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 1, 0 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 0, 1 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 0, 1, 1 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 0, 1 });
    mrc_domain_add_patch(mhd->domain, 1, (int [3]) { 1, 1, 1 });
  } else if (mhd->amr == 11) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 0; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 12) { // 3D
    for (int k = 0; k < 8; k++) {
      for (int j = 0; j < 8; j++) {
	for (int i = 0; i < 8; i++) {
	  mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 13) { // 3D
    /* for (int k = 0; k < 4; k++) { */
    /*   for (int j = 0; j < 4; j++) { */
    /* 	for (int i = 2; i < 4; i++) { */
    /* 	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k }); */
    /* 	} */
    /*   } */
    /* } */
    for (int k = 0; k < 8; k++) {
      for (int j = 0; j < 8; j++) {
	for (int i = 0; i < 1; i++) {
	  mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	}
      }
    }
    for (int k = 0; k < 16; k++) {
      for (int j = 0; j < 16; j++) {
	for (int i = 2; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 4, (int [3]) { i, j, k });
	}
      }
    }
    for (int k = 0; k < 8; k++) {
      for (int j = 0; j < 8; j++) {
	for (int i = 2; i < 8; i++) {
	  mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 14) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 1; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 1; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 6; k++) {
      for (int j = 2; j < 6; j++) {
	    for (int i = 0; i < 4; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 3; k < 4; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 4; k++) {
      for (int j = 3; j < 4; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    //
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 2; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 15) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 1; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 1; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 6; k++) {
      for (int j = 2; j < 6; j++) {
	    for (int i = 0; i < 1; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }

    //
    for (int k = 2; k < 3; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 3; k++) {
      for (int j = 3; j < 5; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 3; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 3; k < 5; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 6; k < 10; k++) {
      for (int j = 6; j < 10; j++) {
	    for (int i = 2; i < 6; i++) {
	      mrc_domain_add_patch(mhd->domain, 4, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 3; k < 5; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 5; k < 6; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 5; k < 6; k++) {
      for (int j = 3; j < 5; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 5; k < 6; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 2; k < 6; k++) {
      for (int j = 2; j < 6; j++) {
	    for (int i = 3; i < 4; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 3; k < 4; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 4; k++) {
      for (int j = 3; j < 4; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    //
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 2; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 16) { // 3D
    for (int k = 0; k < 8; k++) {
      for (int j = 0; j < 8; j++) {
	for (int i = 0; i < 2; i++) {
	  mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	}
      }
    }
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 1; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 17) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 1; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 1; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 6; k++) {
      for (int j = 2; j < 6; j++) {
	    for (int i = 0; i < 1; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }

    //
    for (int k = 2; k < 3; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 3; k++) {
      for (int j = 3; j < 5; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 2; k < 3; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 3; k < 5; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 6; k < 10; k++) {
      for (int j = 6; j < 10; j++) {
	    for (int i = 2; i < 6; i++) {
	      mrc_domain_add_patch(mhd->domain, 4, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 3; k < 5; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 5; k < 6; k++) {
      for (int j = 2; j < 3; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 5; k < 6; k++) {
      for (int j = 3; j < 5; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 5; k < 6; k++) {
      for (int j = 5; j < 6; j++) {
	    for (int i = 1; i < 3; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    //
    for (int k = 2; k < 6; k++) {
      for (int j = 2; j < 6; j++) {
	    for (int i = 3; i < 4; i++) {
	      mrc_domain_add_patch(mhd->domain, 3, (int [3]) { i, j, k });
	    }
      }
    }
    for (int k = 3; k < 4; k++) {
      for (int j = 1; j < 3; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    for (int k = 0; k < 4; k++) {
      for (int j = 3; j < 4; j++) {
	    for (int i = 0; i < 2; i++) {
	      mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	    }
      }
    }

    //
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 2; i < 3; i++) {
	  if ((j == 1 || j == 2) && (k == 1 || k == 2)) {
	    for (int kk = 0; kk < 2; kk++) {
	      for (int jj = 0; jj < 2; jj++) {
		for (int ii = 0; ii < 2; ii++) {
		  mrc_domain_add_patch(mhd->domain, 3,
				       (int [3]) { 2*i+ii, 2*j+jj, 2*k+kk });
		}
	      }
	    }
	  } else {
	    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	  }
	}
      }
    }
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
	for (int i = 3; i < 4; i++) {
	  mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 18) { // 3D
    for (int k = 0; k < 2; k++) {
      for (int i = 0; i < 2; i++) {
	for (int j = 0; j < 2; j++) {
	  mrc_domain_add_patch(mhd->domain, 1, (int [3]) { i, j, k });
	}
      }
    }
  } else if (mhd->amr == 19) { // 3D
    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 4; j++) {
	  if ((i == 1 || i == 2) && (j == 1 || j == 2) && (k == 1 || k == 2)) {
	    for (int kk = 0; kk < 2; kk++) {
	      for (int jj = 0; jj < 2; jj++) {
		for (int ii = 0; ii < 2; ii++) {
		  mrc_domain_add_patch(mhd->domain, 3,
				       (int [3]) { 2*i+ii, 2*j+jj, 2*k+kk });
		}
	      }
	    }
	  } else {
	    mrc_domain_add_patch(mhd->domain, 2, (int [3]) { i, j, k });
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

