
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
  // FIXME, leaked

  struct cuda_mfields *cmflds = cuda_mfields_create();
  cuda_mfields_ctor(cmflds, json);

  mrc_json_t json_cmflds = cuda_mfields_to_json(cmflds); // FIXME, leaked
  int n_patches = mrc_json_get_object_entry_integer(json_cmflds, "n_patches");
  int n_fields = mrc_json_get_object_entry_integer(json_cmflds, "n_fields");

  for (int p = 0; p < n_patches; p++) {
    for (int m = 0; m < n_fields; m++) {
      cuda_mfields_zero_comp_yz(cmflds, m, p);
    }
  }

  cuda_mfields_dump(cmflds, "cmflds.json");

  cuda_mfields_dtor(cmflds);
  cuda_mfields_destroy(cmflds);
}
