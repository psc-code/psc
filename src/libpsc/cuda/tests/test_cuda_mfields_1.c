
#include "cuda_iface.h"

#include "json-builder.h"

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
// main

int
main(void)
{
  mrc_json_t json = mrc_json_parse("{                                           "
				   "  \"info\" : {                              "
				   "    \"n_patches\" : 1,                      "
				   "    \"n_fields\" : 6,                       "
				   "    \"ldims\" : [ 1, 8, 8 ],                "
				   "    \"ib\" : [ 0, -2, -2 ],                 "
				   "    \"im\" : [ 1, 12, 12 ],                 "
				   "    \"dx\" : [ 1.0, 10.0, 10.0 ],           "
				   "  }                                         "
				   "}                                           ");
  mrc_json_print(json, 0);

  struct cuda_mfields *cmflds = cuda_mfields_create();
  cuda_mfields_ctor(cmflds, json);
  cuda_mfields_dump(cmflds);

#if 0
  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  unsigned int n_prts_by_patch[n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, n_patches, n_prts_by_patch, json_info);
  printf("added particles\n");
  cuda_mparticles_dump_by_patch(cmprts, n_prts_by_patch);

  cuda_mparticles_setup_internals(cmprts);
  printf("set up internals\n");
  cuda_mparticles_dump(cmprts);

  printf("get_particles_test\n");
  get_particles_test(cmprts, n_patches, n_prts_by_patch);
#endif
  
  cuda_mfields_dtor(cmflds);
  cuda_mfields_destroy(cmflds);
}
