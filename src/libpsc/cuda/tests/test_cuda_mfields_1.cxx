
#include "cuda_mfields.h"

#include "fields.hxx"

#include <cmath>

enum { // FIXME, duplicated
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
  NR_FIELDS,
};

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mrc_profile.h>

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

// ----------------------------------------------------------------------
// init_wave

double
init_wave(double x, double y, int m)
{
  double kx = 2. * M_PI / 80., ky = 2. * M_PI / 80.;
  switch (m) {
  case EX: return   1./sqrtf(2.) * sin(kx * x + ky * y);
  case EY: return - 1./sqrtf(2.) * sin(kx * x + ky * y);
  case HZ: return sin(kx * x + ky * y);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// main

int
main(void)
{
  using fields_t = fields_single_t;
  using real_t = fields_t::real_t;
  
  mrc_json_t json = mrc_json_parse("{                                           "
				   "  \"info\" : {                              "
				   "    \"n_fields\" : 9,                       "
				   "    \"ib\" : [ 0, -2, -2 ],                 "
				   "    \"im\" : [ 1, 12, 12 ],                 "
				   "    \"dx\" : [ 1.0, 10.0, 10.0 ],           "
				   "  }                                         "
				   "}                                           ");
  mrc_json_print(json, 0);
  // FIXME, leaked
  Grid_t grid{};
  grid.ldims = { 1, 8, 8 };

  Grid_t::Patch patch{};
  patch.xb = { 0., 0., 0. };
  grid.patches.push_back(patch);

  struct cuda_mfields *cmflds = new cuda_mfields(grid, json);

  mrc_json_t json_cmflds = cuda_mfields_to_json(cmflds); // FIXME, leaked
  int n_patches = cmflds->n_patches;
  int n_fields = mrc_json_get_object_entry_integer(json_cmflds, "n_fields");
  double dx[3]; mrc_json_get_object_entry_double3(json_cmflds, "dx", dx);

  for (int m = 0; m < n_fields; m++) {
    cuda_mfields_zero_comp_yz(cmflds, m);
  }

  cuda_mfields_dump(cmflds, "cmflds.json");

  fields_t flds = cuda_mfields_get_host_fields(cmflds);
  Fields3d<fields_t> F(flds);
  for (int p = 0; p < n_patches; p++) {
    for (int k = flds.ib[2]; k < flds.ib[2] + flds.im[2]; k++) {
      for (int j = flds.ib[1]; j < flds.ib[1] + flds.im[1]; j++) {
	for (int i = flds.ib[0]; i < flds.ib[0] + flds.im[0]; i++) {
	  real_t x_cc = (i + .5) * dx[0];
	  real_t y_cc = (j + .5) * dx[1];
	  real_t x_nc = i * dx[0];
	  real_t y_nc = j * dx[1];
	  F(EX, i,j,k) = init_wave(x_cc, y_nc, EX);
	  F(EY, i,j,k) = init_wave(x_nc, y_cc, EY);
	  F(EZ, i,j,k) = init_wave(x_nc, y_nc, EZ);
	  F(HX, i,j,k) = init_wave(x_nc, y_cc, HX);
	  F(HY, i,j,k) = init_wave(x_cc, y_nc, HY);
	  F(HZ, i,j,k) = init_wave(x_cc, y_cc, HZ);
	}
      }
    }
    cuda_mfields_copy_to_device(cmflds, p, flds, 0, n_fields);
  }
  flds.dtor();

  cuda_mfields_dump(cmflds, "cmflds_wave.json");

  float dt = dx[1];
  cuda_push_fields_E_yz(cmflds, .5 * dt);
  cuda_push_fields_H_yz(cmflds, dt);
  cuda_push_fields_E_yz(cmflds, .5 * dt);
  
  cuda_mfields_dump(cmflds, "cmflds_wave_1.json");
  
  delete cmflds;
}
