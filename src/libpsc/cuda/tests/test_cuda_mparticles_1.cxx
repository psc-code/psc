
#include "cuda_iface.h"

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
  mrc_json_t json_info = *(mrc_json_t *) ctx;

  int ldims[3];
  double dx[3];
  mrc_json_get_object_entry_int3(json_info, "ldims", ldims);
  mrc_json_get_object_entry_double3(json_info, "dx", dx);

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
				     int n_patches,
				     unsigned int *n_prts_by_patch,
				     mrc_json_t json_info)
{
  int ldims[3];
  mrc_json_get_object_entry_int3(json_info, "ldims", ldims);

  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = ldims[0] * ldims[1] * ldims[2];
  }

  cmprts->reserve_all(n_prts_by_patch);
  cmprts->resize_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    cmprts->set_particles(n_prts_by_patch[p], off, set_particle_test_1,
			  &json_info);
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
get_particles_test(struct cuda_mparticles *cmprts, int n_patches,
		   unsigned int *n_prts_by_patch)
{
  unsigned int off = 0;
  for (int p = 0; p < n_patches; p++) {
    cmprts->get_particles(n_prts_by_patch[p], off, get_particle, NULL);
    off += n_prts_by_patch[p];
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
				   "    \"ldims\" : [ 1, 4, 2 ],                "
				   "    \"bs\" : [ 1, 1, 1 ],                   "
				   "    \"dx\" : [ 1.0, 10.0, 10.0 ],           "
				   "    \"xb_by_patch\" : [ [ 0.0, 0.0, 0.0 ] ],"
				   "    \"fnqs\" : 1.0,                         "
				   "    \"eta\" : 1.0,                          "
				   "    \"dt\" : 1.0,                           "
				   "    \"kind_q\" : [ -1.0, 1.0 ],             "
				   "    \"kind_m\" : [  1.0, 25.0 ]             "
				   "  }                                         "
				   "}                                           ");
  mrc_json_print(json, 0);

  Grid<double> grid;
  struct cuda_mparticles *cmprts = new cuda_mparticles(grid, json);

  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  int n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  unsigned int n_prts_by_patch[n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, n_patches, n_prts_by_patch, json_info);
  printf("added particles\n");
  cmprts->dump_by_patch(n_prts_by_patch);

  cmprts->setup_internals();
  printf("set up internals\n");
  cmprts->dump();

  printf("get_particles_test\n");
  get_particles_test(cmprts, n_patches, n_prts_by_patch);

  delete cmprts;
}
