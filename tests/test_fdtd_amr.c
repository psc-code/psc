
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_ddc.h>
#include <mrc_io.h>
#include <mrctest.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define ARRAY_SIZE(a) (sizeof(a) / sizeof(a[0]))

// ----------------------------------------------------------------------
// x-y plane, Xx: interior point, Oo: ghost point

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

// Ex
// +-o-X-o-+-x-O-x-+
// |       |   |   |
// |       +-x-+-x-+
// |       |   |   |
// +-o-X-o-+-x-O-x-+

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

#define F3 MRC_M3

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

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EX[2] = {
  // FIXME, 3D
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 1, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EY[2] = {
  // FIXME, 3D
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 1, 0, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_EZ[4] = {
  { .dx = { 0, 0, 0 }, .val = .25f },
  { .dx = { 1, 0, 0 }, .val = .25f },
  { .dx = { 0, 1, 0 }, .val = .25f },
  { .dx = { 1, 1, 0 }, .val = .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HX[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 1, 0, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HY[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 1, 0 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil_entry stencil_coarse_HZ[2] = {
  { .dx = { 0, 0, 0 }, .val = .5f },
  { .dx = { 0, 0, 1 }, .val = .5f },
};

static struct mrc_ddc_amr_stencil stencils_coarse[NR_COMPS] = {
  [EX] = { stencil_coarse_EX, ARRAY_SIZE(stencil_coarse_EX) },
  [EY] = { stencil_coarse_EY, ARRAY_SIZE(stencil_coarse_EY) },
  [EZ] = { stencil_coarse_EZ, ARRAY_SIZE(stencil_coarse_EZ) },
  [HX] = { stencil_coarse_HX, ARRAY_SIZE(stencil_coarse_HX) },
  [HY] = { stencil_coarse_HY, ARRAY_SIZE(stencil_coarse_HY) },
  [HZ] = { stencil_coarse_HZ, ARRAY_SIZE(stencil_coarse_HZ) },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EX[6] = {
  // FIXME, 3D
  { .dx = {  0, -1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { +1, -1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = {  0, +1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (1.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EY[6] = {
  // FIXME, 3D
  { .dx = { -1,  0,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0,  0,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1,  0,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = { -1, +1,  0 }, .val = (1.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (1.f/8.f) * 2.f },
  { .dx = { +1, +1,  0 }, .val = (1.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_EZ[9] = {
  // FIXME, 3D
  { .dx = { -1, -1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = {  0, -1,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { +1, -1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = { -1,  0,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f  },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { -1, +1,  0 }, .val = (2.f/8.f) * .25f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * .5f  },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .25f },
};

static struct mrc_ddc_amr_stencil_entry stencil_fine_HX[6] = {
  // FIXME, 3D
  { .dx = { -1,  0,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { -1, +1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .5f },
};
	  
static struct mrc_ddc_amr_stencil_entry stencil_fine_HY[6] = {
  // FIXME, 3D
  { .dx = {  0, -1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { +1, -1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * .5f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * .5f },
};
	  
static struct mrc_ddc_amr_stencil_entry stencil_fine_HZ[4] = {
  // FIXME, 3D
  { .dx = {  0,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1,  0,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = {  0, +1,  0 }, .val = (2.f/8.f) * 1.f },
  { .dx = { +1, +1,  0 }, .val = (2.f/8.f) * 1.f },
};

static struct mrc_ddc_amr_stencil stencils_fine[NR_COMPS] = {
  [EX] = { stencil_fine_EX, ARRAY_SIZE(stencil_fine_EX) },
  [EY] = { stencil_fine_EY, ARRAY_SIZE(stencil_fine_EY) },
  [EZ] = { stencil_fine_EZ, ARRAY_SIZE(stencil_fine_EZ) },
  [HX] = { stencil_fine_HX, ARRAY_SIZE(stencil_fine_HX) },
  [HY] = { stencil_fine_HY, ARRAY_SIZE(stencil_fine_HY) },
  [HZ] = { stencil_fine_HZ, ARRAY_SIZE(stencil_fine_HZ) },
};

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

#define mrc_domain_amr(domain) mrc_to_subobj(domain, struct mrc_domain_amr)

static void
step_fdtd(struct mrc_m3 *fld, struct mrc_ddc *ddc_E, struct mrc_ddc *ddc_H)
{
  struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
  int ldims[3], nr_levels;
  mrc_domain_get_param_int3(fld->domain, "m", ldims);
  mrc_domain_get_nr_levels(fld->domain, &nr_levels);
  float dx = 1. / (ldims[0] << (nr_levels - 1));
  float dt = dx / sqrt(2.);

  mrc_ddc_amr_apply(ddc_H, fld);

  mrc_m3_foreach_patch(fld, p) {
    struct mrc_m3_patch *fldp = mrc_m3_patch_get(fld, p);
    mrc_crds_patch_get(crds, p);
    float dx = MRC_MCRDX(crds, 1) - MRC_MCRDX(crds, 0); // FIXME
    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = .5 * dt / dy;
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
    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = .5 * dt / dy;
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
    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = .5 * dt / dy;
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
    float dy = MRC_MCRDY(crds, 1) - MRC_MCRDY(crds, 0);
    float cnx = .5 * dt / dx;
    float cny = .5 * dt / dy;
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

float
func1(float x, float y, int m)
{
  float kx = 2. * M_PI;
  return sin(.5 + kx * x);
}

float
func2(float x, float y, int m)
{
  float kx = 2. * M_PI, ky = 2. * M_PI;
  return sin(.5 + kx * x) * cos(.5 + ky * y);
}

float
func3(float x, float y, int m)
{
  float kx = 2. * M_PI, ky = 2. * M_PI;
  switch (m) {
  case EX: return   1./sqrtf(2.) * sin(kx * x + ky * y);
  case EY: return - 1./sqrtf(2.) * sin(kx * x + ky * y);
  case HZ: return sin(kx * x + ky * y);
  default: return 0.;
  }
}

float (*func)(float, float, int) = func3;

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
  mrctest_set_amr_domain_4(domain);

  mrc_domain_setup(domain);
  mrc_domain_plot(domain);

  // create and fill a field

  struct mrc_m3 *fld = mrc_domain_m3_create(domain);
  mrc_m3_set_name(fld, "fld");
  mrc_m3_set_nr_comps(fld, NR_COMPS);
  mrc_m3_set_sw(fld, 3);
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
      MRC_M3(fldp, EX, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, EY, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, EZ, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, HX, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, HY, ix,iy,iz) = 1.f / 0.f;
      MRC_M3(fldp, HZ, ix,iy,iz) = 1.f / 0.f;
    } mrc_m3_foreach_end;
#endif
    mrc_m3_foreach(fldp, ix,iy,iz, 0, 1) {
      float x_cc = MRC_MCRDX(crds, ix);
      float y_cc = MRC_MCRDY(crds, iy);
      float x_nc = .5f * (MRC_MCRDX(crds, ix-1) + MRC_MCRDX(crds, ix));
      float y_nc = .5f * (MRC_MCRDY(crds, iy-1) + MRC_MCRDY(crds, iy));
      if (!mrc_domain_is_ghost(domain, (int[]) { 0, 1, 1 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, EX, ix,iy,iz) = func(x_cc, y_nc, EX);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 1, 0, 1 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, EY, ix,iy,iz) = func(x_nc, y_cc, EY);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 1, 1, 0 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, EZ, ix,iy,iz) = func(x_nc, y_nc, EZ);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 1, 0, 0 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, HX, ix,iy,iz) = func(x_nc, y_cc, HX);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 0, 1, 0 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, HY, ix,iy,iz) = func(x_cc, y_nc, HY);
      }
      if (!mrc_domain_is_ghost(domain, (int[]) { 0, 0, 1 }, p, (int[]) { ix, iy, iz })) {
	MRC_M3(fldp, HZ, ix,iy,iz) = func(x_cc, y_cc, HZ);
      }
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(fld);
    mrc_crds_patch_put(crds);
  }

  struct mrc_ddc *ddc_E = mrc_ddc_create(mrc_domain_comm(domain));
  mrc_ddc_set_type(ddc_E, "amr");
  mrc_ddc_set_domain(ddc_E, domain);
  mrc_ddc_set_param_int(ddc_E, "size_of_type", sizeof(float));
  mrc_ddc_set_param_int(ddc_E, "sw", fld->_sw.vals[0]); // FIXME, 1,2 maybe different
  mrc_ddc_setup(ddc_E);
  mrc_ddc_amr_set_by_stencil(ddc_E, EX, 2, (int[]) { 0, 1, 1 }, &stencils_coarse[EX], &stencils_fine[EX]);
  mrc_ddc_amr_set_by_stencil(ddc_E, EY, 2, (int[]) { 1, 0, 1 }, &stencils_coarse[EY], &stencils_fine[EY]);
  mrc_ddc_amr_set_by_stencil(ddc_E, EZ, 2, (int[]) { 1, 1, 0 }, &stencils_coarse[EZ], &stencils_fine[EZ]);
  mrc_ddc_amr_assemble(ddc_E);

  struct mrc_ddc *ddc_H = mrc_ddc_create(mrc_domain_comm(domain));
  mrc_ddc_set_type(ddc_H, "amr");
  mrc_ddc_set_domain(ddc_H, domain);
  mrc_ddc_set_param_int(ddc_H, "size_of_type", sizeof(float));
  mrc_ddc_set_param_int(ddc_H, "sw", fld->_sw.vals[0]); // FIXME
  mrc_ddc_setup(ddc_H);
  mrc_ddc_amr_set_by_stencil(ddc_H, HX, 2, (int[]) { 1, 0, 0 }, &stencils_coarse[HX], &stencils_fine[HX]);
  mrc_ddc_amr_set_by_stencil(ddc_H, HY, 2, (int[]) { 0, 1, 0 }, &stencils_coarse[HY], &stencils_fine[HY]);
  mrc_ddc_amr_set_by_stencil(ddc_H, HZ, 2, (int[]) { 0, 0, 1 }, &stencils_coarse[HZ], &stencils_fine[HZ]);
  mrc_ddc_amr_assemble(ddc_H);

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

  for (int n = 0; n <= 100; n++) {
    mrc_io_open(io, "w", n+2, n+2);
    mrc_m3_write(fld, io);
    mrc_io_close(io);

    step_fdtd(fld, ddc_E, ddc_H);
  }

  mrc_io_destroy(io);

  mrc_ddc_destroy(ddc_E);
  mrc_ddc_destroy(ddc_H);

  mrc_m3_destroy(fld);

  mrc_domain_destroy(domain);

  MPI_Finalize();
}
