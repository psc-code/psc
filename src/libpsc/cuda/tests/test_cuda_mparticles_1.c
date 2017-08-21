
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
// cuda_domain_info_ctor_test_1

void
cuda_domain_info_ctor_test_1(struct cuda_domain_info *info)
{
  info->n_patches = 1;
  info->ldims[0] = 1; info->ldims[1] = 4; info->ldims[2] = 2;
  info->bs[0] = 1; info->bs[1] = 1; info->bs[2] = 1;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;

  info->xb_by_patch = calloc(info->n_patches, sizeof(*info->xb_by_patch));
  for (int p = 0; p < info->n_patches; p++) {
    info->xb_by_patch[p][0] = 0.;
    info->xb_by_patch[p][1] = 0.;
    info->xb_by_patch[p][2] = 0.;
  }
};

// ----------------------------------------------------------------------
// cuda_domain_info_dtor

void
cuda_domain_info_dtor(struct cuda_domain_info *info)
{
  free(info->xb_by_patch);
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
  struct cuda_domain_info *info = ctx;

  int k = n % info->ldims[2];
  n /= info->ldims[2];
  int j = n % info->ldims[1];
  n /= info->ldims[1];
  int i = n;

  prt->xi[0] = info->dx[0] * (i + .5f);
  prt->xi[1] = info->dx[1] * (j + .5f);
  prt->xi[2] = info->dx[2] * (k + .5f);
  prt->pxi[0] = i;
  prt->pxi[1] = j;
  prt->pxi[2] = k;
  prt->kind = 0;
  prt->qni_wni = 1.;
}

void
cuda_mparticles_add_particles_test_1(struct cuda_mparticles *cmprts,
				     struct cuda_domain_info *info,
				     unsigned int *n_prts_by_patch)
{
  for (int p = 0; p < info->n_patches; p++) {
    n_prts_by_patch[p] = info->ldims[0] * info->ldims[1] * info->ldims[2];
  }

  cuda_mparticles_reserve_all(cmprts, n_prts_by_patch);
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < info->n_patches; p++) {
    cuda_mparticles_set_particles(cmprts, n_prts_by_patch[p], off,
				  set_particle_test_1, info);
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
		   struct cuda_domain_info *info,
		   unsigned int *n_prts_by_patch)
{
  unsigned int off = 0;
  for (int p = 0; p < info->n_patches; p++) {
    cuda_mparticles_get_particles(cmprts, n_prts_by_patch[p], off,
				  get_particle, info);
    off += n_prts_by_patch[p];
  }
}

// ----------------------------------------------------------------------
// main

int
main(void)
{
  struct cuda_mparticles *cmprts = cuda_mparticles_create();

  struct cuda_domain_info info;
  cuda_domain_info_ctor_test_1(&info);
  cuda_mparticles_set_domain_info(cmprts, &info);
  cuda_domain_info_dtor(&info);

  cuda_mparticles_setup(cmprts);

  unsigned int n_prts_by_patch[info.n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, &info, n_prts_by_patch);
  printf("added particles\n");
  cuda_mparticles_dump_by_patch(cmprts, n_prts_by_patch);

  cuda_mparticles_setup_internals(cmprts);
  printf("set up internals\n");
  cuda_mparticles_dump(cmprts);

  printf("get_particles_test\n");
  get_particles_test(cmprts, &info, n_prts_by_patch);

  cuda_mparticles_destroy(cmprts);
}
