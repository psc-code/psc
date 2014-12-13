
#include "ggcm_mhd_step.h" // FIXME

// FIXME, storing these ddc's in ggcm_mhd is kinda ugly, and they
// still don't get cleaned up, anyway

// ----------------------------------------------------------------------
// correct_E

// FIXME -> header
#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

void mrc_domain_get_neighbor_patch_same(struct mrc_domain *domain, int p,
					int dx[3], int *p_nei);
void mrc_domain_get_neighbor_patch_fine(struct mrc_domain *domain, int gp,
					int dir[3], int off[3], int *gp_nei);

void mrc_domain_find_valid_point_fine(struct mrc_domain *domain, int ext[3], int gp, int i[3],
				      int *gp_nei, int j[3]);
void mrc_ddc_amr_add_diagonal_one(struct mrc_ddc *ddc, int gp, int m, int i[3]);

// ----------------------------------------------------------------------
// create_ddc_amr_E

static struct mrc_ddc *
create_ddc_amr_E(struct ggcm_mhd *mhd)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  mrc_ddc_set_param_int(ddc, "n_comp", 3);
  mrc_ddc_setup(ddc);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int sw[3] = {};
  int ext[3] = { 1, 1, 0 }; // EZ

  int nr_patches;
  mrc_domain_get_patches(mhd->domain, &nr_patches);
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mhd->domain, p, &info);
    int gp = info.global_patch, *ldims = info.ldims;

    int i[3];
    for (i[2] = -sw[2]; i[2] < ldims[2] + ext[2] + sw[2]; i[2]++) {
      for (i[1] = -sw[1]; i[1] < ldims[1] + ext[1] + sw[1]; i[1]++) {
	for (i[0] = -sw[0]; i[0] < ldims[0] + ext[0] + sw[0]; i[0]++) {
	  // 2D EZ only
	  if (i[0] > 0 && i[0] < ldims[0] &&
	      i[1] > 0 && i[1] < ldims[1]) {
	    // truly interior point "x", ie., not on boundary
	    // X---X---X---X 
	    // |   |   |   |
	    // X---x---x---X
	    // |   |   |   |
	    // X---x---x---X
	    // |   |   |   |
	    // X---X---X---X 
	    mrc_ddc_amr_add_diagonal_one(ddc, gp, 2, i);
	    continue;
	  }

	  // now we're only looking at EZ field values that are on edges on the boundary of
	  // this patch

	  // FIXME, if there's another patch that really owns the value, it should be
	  // overwritten from there (though it gets calculated redundantly, so that's why
	  // we get away without doing it)

	  // If we're bordering a fine patch so that this EZ value is also on the fine grid,
	  // replace this coarse grid value by the one from the fine grid in the same place

	  int gp_nei, j[3];
	  mrc_domain_find_valid_point_fine(mhd->domain, ext, gp,
					   (int[]) { 2*i[0], 2*i[1], 2*i[2] }, &gp_nei, j);
	  if (gp_nei >= 0) {
	    mrc_ddc_amr_add_value(ddc, gp, 2, i, gp_nei, 2, j, 1.f);
	    continue;
	  }

	  // FIXME not handled -- should be fixed, though (see FIXME above) not strictly necessary
	  // MHERE;
	}
      }
    }
  }

  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// ----------------------------------------------------------------------
// correct_E

static void
correct_E(struct ggcm_mhd *mhd, struct mrc_fld *E, int l, int r)
{
  if (!mhd->ddc_amr_E) {
    // FIXME, leaked
    mhd->ddc_amr_E = create_ddc_amr_E(mhd);
  }

  mrc_ddc_amr_apply(mhd->ddc_amr_E, E);
}

// ======================================================================
// correct fluxes
// ======================================================================

static struct mrc_ddc_amr_stencil_entry stencil_fine_flux_x[] = {
  // FIXME, 3D
  { .dx = {  0,  0,  0 }, .val = .5f },
  { .dx = {  0, +1,  0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil stencils_fine_flux_x = {
  stencil_fine_flux_x, ARRAY_SIZE(stencil_fine_flux_x)
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_flux_y[] = {
  // FIXME, 3D
  { .dx = {  0,  0,  0 }, .val = .5f },
  { .dx = { +1,  0,  0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil stencils_fine_flux_y = {
  stencil_fine_flux_y, ARRAY_SIZE(stencil_fine_flux_y)
};

// ----------------------------------------------------------------------
// create_ddc_amr_flux_x

static struct mrc_ddc *
create_ddc_amr_flux_x(struct ggcm_mhd *mhd)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  // FIXME!!!
  if (strcmp(ggcm_mhd_step_type(mhd->step), "vl") == 0) {
    mrc_ddc_set_param_int(ddc, "n_comp", 5);
  } else if (strcmp(ggcm_mhd_step_type(mhd->step), "vlct") == 0) {
    mrc_ddc_set_param_int(ddc, "n_comp", 8);
  } else {
    assert(0);
  }
  mrc_ddc_setup(ddc);
  for (int m = 0; m < 5; m++) {
    mrc_ddc_amr_set_by_stencil(ddc, m, 2, (int[]) { 1, 0, 0 },
			       NULL, &stencils_fine_flux_x);
  }
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// ----------------------------------------------------------------------
// create_ddc_amr_flux_y

static struct mrc_ddc *
create_ddc_amr_flux_y(struct ggcm_mhd *mhd)
{
  struct mrc_ddc *ddc = mrc_ddc_create(mrc_domain_comm(mhd->domain));
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, mhd->domain);
  mrc_ddc_set_param_int(ddc, "size_of_type", mhd->fld->_size_of_type);
  mrc_ddc_set_param_int3(ddc, "sw", mrc_fld_spatial_sw(mhd->fld));
  if (strcmp(ggcm_mhd_step_type(mhd->step), "vl") == 0) {
    mrc_ddc_set_param_int(ddc, "n_comp", 5);
  } else if (strcmp(ggcm_mhd_step_type(mhd->step), "vlct") == 0) {
    mrc_ddc_set_param_int(ddc, "n_comp", 8);
  } else {
    assert(0);
  }
  mrc_ddc_setup(ddc);
  for (int m = 0; m < 5; m++) {
    mrc_ddc_amr_set_by_stencil(ddc, m, 2, (int[]) { 0, 1, 0 },
			       NULL, &stencils_fine_flux_y);
  }
  mrc_ddc_amr_assemble(ddc);

  return ddc;
}

// ----------------------------------------------------------------------
// correct_fluxes

static void
correct_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *fluxes[3])
{
  if (!mhd->ddc_amr_flux_x) {
    // FIXME, leaked
    mhd->ddc_amr_flux_x = create_ddc_amr_flux_x(mhd);
  }
  if (!mhd->ddc_amr_flux_y) {
    // FIXME, leaked
    mhd->ddc_amr_flux_y = create_ddc_amr_flux_y(mhd);
  }

  mrc_ddc_amr_apply(mhd->ddc_amr_flux_x, fluxes[0]);
  mrc_ddc_amr_apply(mhd->ddc_amr_flux_y, fluxes[1]);
}

