
#include "cuda_iface.h"

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
  mrc_json_t json = mrc_json_parse("{                                           "
				   "  \"info\" : {                              "
				   "    \"n_patches\" : 1,                      "
				   "    \"n_fields\" : 9,                       "
				   "    \"ldims\" : [ 1, 8, 8 ],                "
				   "    \"ib\" : [ 0, -2, -2 ],                 "
				   "    \"im\" : [ 1, 12, 12 ],                 "
				   "    \"dx\" : [ 1.0, 10.0, 10.0 ],           "
				   "  }                                         "
				   "}                                           ");
  mrc_json_print(json, 0);
  // FIXME, leaked
  struct cuda_mfields *cmflds = cuda_mfields_create();
  cuda_mfields_ctor(cmflds, json);

  mrc_json_t json_cmflds = cuda_mfields_to_json(cmflds); // FIXME, leaked
  int n_patches = mrc_json_get_object_entry_integer(json_cmflds, "n_patches");
  int n_fields = mrc_json_get_object_entry_integer(json_cmflds, "n_fields");
  double dx[3]; mrc_json_get_object_entry_double3(json_cmflds, "dx", dx);

  for (int m = 0; m < n_fields; m++) {
    cuda_mfields_zero_comp_yz(cmflds, m);
  }

  cuda_mfields_dump(cmflds, "cmflds.json");

  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  for (int p = 0; p < n_patches; p++) {
    for (int k = flds.ib[2]; k < flds.ib[2] + flds.im[2]; k++) {
      for (int j = flds.ib[1]; j < flds.ib[1] + flds.im[1]; j++) {
	for (int i = flds.ib[0]; i < flds.ib[0] + flds.im[0]; i++) {
	  fields_single_real_t x_cc = (i + .5) * dx[0];
	  fields_single_real_t y_cc = (j + .5) * dx[1];
	  fields_single_real_t x_nc = i * dx[0];
	  fields_single_real_t y_nc = j * dx[1];
	  _F3_S(flds, EX, i,j,k) = init_wave(x_cc, y_nc, EX);
	  _F3_S(flds, EY, i,j,k) = init_wave(x_nc, y_cc, EY);
	  _F3_S(flds, EZ, i,j,k) = init_wave(x_nc, y_nc, EZ);
	  _F3_S(flds, HX, i,j,k) = init_wave(x_nc, y_cc, HX);
	  _F3_S(flds, HY, i,j,k) = init_wave(x_cc, y_nc, HY);
	  _F3_S(flds, HZ, i,j,k) = init_wave(x_cc, y_cc, HZ);
	}
      }
    }
    cuda_mfields_copy_to_device(cmflds, p, flds, 0, n_fields);
  }
  fields_single_t_dtor(&flds);

  cuda_mfields_dump(cmflds, "cmflds_wave.json");

  float dt = dx[1];
  cuda_push_fields_E_yz(cmflds, .5 * dt);
  cuda_push_fields_H_yz(cmflds, dt);
  cuda_push_fields_E_yz(cmflds, .5 * dt);
  
  cuda_mfields_dump(cmflds, "cmflds_wave_1.json");
  
  cuda_mfields_dtor(cmflds);
  cuda_mfields_destroy(cmflds);
}
