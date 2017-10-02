
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
// cuda_mparticles_add_particles_test_1
//
// add 1 particle at the center of each cell, in the "wrong" order in each
// patch (C order, but to get them ordered by block, they need to be reordered
// into Fortran order, a.k.a., this will exercise the initial sorting

static void
set_particle_test_1(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  mrc_json_t json = *(mrc_json_t *) ctx;
  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int ldims[3];
  double dx[3];

  mrc_json_t json_ldims = mrc_json_get_object_entry(json_info, "ldims");
  mrc_json_t json_dx = mrc_json_get_object_entry(json_info, "dx");
  for (int d = 0; d < 3; d++) {
    ldims[d] = mrc_json_get_array_entry_integer(json_ldims, d);
    dx[d] = mrc_json_get_array_entry_double(json_dx, d);
  }

  int k = n % ldims[2];
  n /= ldims[2];
  int j = n % ldims[1];
  n /= ldims[1];
  int i = n;

  prt->xi[0] = dx[0] * (i + .5f);
  prt->xi[1] = dx[1] * (j + .5f);
  prt->xi[2] = dx[2] * (k + .5f);
  prt->pxi[0] = i;
  prt->pxi[1] = j;
  prt->pxi[2] = k;
  prt->kind = 0;
  prt->qni_wni = 1.;
}

void
cuda_mparticles_add_particles_test_1(struct cuda_mparticles *cmprts,
				     mrc_json_t json,
				     unsigned int *n_prts_by_patch)
{
  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  int ldims[3];

  mrc_json_t json_ldims = mrc_json_get_object_entry(json_info, "ldims");
  for (int d = 0; d < 3; d++) {
    ldims[d] = mrc_json_get_array_entry_integer(json_ldims, d);
  }
  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = ldims[0] * ldims[1] * ldims[2];
  }

  cuda_mparticles_reserve_all(cmprts, n_prts_by_patch);
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    cuda_mparticles_set_particles(cmprts, n_prts_by_patch[p], off,
				  set_particle_test_1, &json);
    off += n_prts_by_patch[p];
  }
}

// ----------------------------------------------------------------------
// get_particles_test

static void
get_particle(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  printf("prt[%d] xi %g %g %g // pxi %g %g %g // kind %d // qni_wni %g\n",
	 n, prt->xi[0], prt->xi[1], prt->xi[2],
	 prt->pxi[0], prt->pxi[1], prt->pxi[2],
	 prt->kind, prt->qni_wni);
}

void
get_particles_test(struct cuda_mparticles *cmprts,
		   mrc_json_t json,
		   unsigned int *n_prts_by_patch)
{
  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    cuda_mparticles_get_particles(cmprts, n_prts_by_patch[p], off,
				  get_particle, NULL);
    off += n_prts_by_patch[p];
  }
}

// ----------------------------------------------------------------------
// main

int
main(void)
{
  struct cuda_mparticles *cmprts = cuda_mparticles_create();

  mrc_json_t json = mrc_json_parse("{                                           "
				   "  \"info\" : {                              "
				   "    \"n_patches\" : 1,                      "
				   "    \"ldims\" : [ 1, 4, 2 ],                "
				   "    \"bs\" : [ 1, 1, 1 ],                   "
				   "    \"dx\" : [ 1.0, 10.0, 10.0 ],           "
				   "    \"xb_by_patch\" : [ [ 0.0, 0.0, 0.0 ] ] "
				   "  }                                         "
				   "}                                           ");
  mrc_json_print(json, 0);

  cuda_mparticles_set_domain_info(cmprts, json);

  cuda_mparticles_setup(cmprts);

  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  unsigned int n_prts_by_patch[n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, json, n_prts_by_patch);
  printf("added particles\n");
  cuda_mparticles_dump_by_patch(cmprts, n_prts_by_patch);

  cuda_mparticles_setup_internals(cmprts);
  printf("set up internals\n");
  cuda_mparticles_dump(cmprts);

  printf("get_particles_test\n");
  get_particles_test(cmprts, json, n_prts_by_patch);

  cuda_mparticles_destroy(cmprts);
}
