
#include "cuda_iface.h"

#include <stdio.h>
#include <assert.h>

// ----------------------------------------------------------------------
// cuda_domain_info_set_test_1

void
cuda_domain_info_set_test_1(struct cuda_domain_info *info)
{
  info->n_patches = 1;
  info->ldims[0] = 1; info->ldims[1] = 4; info->ldims[2] = 2;
  info->bs[0] = 1; info->bs[1] = 1; info->bs[2] = 1;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;
};

// ----------------------------------------------------------------------
// cuda_mparticles_add_particles_test_1
//
// add 1 particle at the center of each cell, in the "wrong" order in each
// patch (C order, but to get them ordered by block, they need to be reordered
// into Fortran order, a.k.a., this will exercise the initial sorting

static void
get_particle_test_1(struct cuda_mparticles_prt *prt, int n, void *ctx)
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

  cuda_mparticles_alloc(cmprts, n_prts_by_patch);

  unsigned int off = 0;
  for (int p = 0; p < info->n_patches; p++) {
    cuda_mparticles_set_particles(cmprts, n_prts_by_patch[p], off,
				  get_particle_test_1, info);
    off += n_prts_by_patch[p];
  }
}

#if 0
static int
get_block_idx(struct cuda_mparticles *cmprts, int n, int p)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  float *b_dxi = cmprts->b_dxi;
  int *b_mx = cmprts->b_mx;
  
  float4 xi4 = d_xi4[n];
  unsigned int block_pos_y = (int) floor(xi4.y * b_dxi[1]);
  unsigned int block_pos_z = (int) floor(xi4.z * b_dxi[2]);

  int bidx;
  if (block_pos_y >= b_mx[1] || block_pos_z >= b_mx[2]) {
    bidx = -1;
  } else {
    bidx = (p * b_mx[2] + block_pos_z) * b_mx[1] + block_pos_y;
  }

  return bidx;
}

static void
cuda_mparticles_check_in_patch_unordered(struct cuda_mparticles *cmprts,
					 unsigned int *nr_prts_by_patch)
{
  unsigned int off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    for (int n = 0; n < nr_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(cmprts, off + n, p);
      assert(bidx >= 0 && bidx <= cmprts->n_blocks);
    }
    off += nr_prts_by_patch[p];
  }

  assert(off == cmprts->n_prts);
}

static void
cuda_mparticles_check_bidx_id_unordered(struct cuda_mparticles *cmprts,
					unsigned int *n_prts_by_patch)
{
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);

  unsigned int off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(cmprts, off + n, p);
      assert(bidx == d_bidx[off+n]);
      assert(off+n == d_id[off+n]);
    }
    off += n_prts_by_patch[p];
  }

  assert(off == cmprts->n_prts);
}
#endif

int
main(void)
{
  struct cuda_mparticles *cmprts = cuda_mparticles_create();

  struct cuda_domain_info info;
  cuda_domain_info_set_test_1(&info);
  cuda_mparticles_set_domain_info(cmprts, &info);

  unsigned int n_prts_by_patch[info.n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, &info, n_prts_by_patch);
  printf("added particles\n");
  cuda_mparticles_dump(cmprts);
#if 0
  cuda_mparticles_check_in_patch_unordered(cmprts, n_prts_by_patch);

  cuda_mparticles_find_block_indices_ids(cmprts, n_prts_by_patch);
  printf("find bidx, id\n");
  cuda_mparticles_dump(cmprts);
  cuda_mparticles_check_bidx_id_unordered(cmprts, n_prts_by_patch);

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::stable_sort_by_key(d_bidx, d_bidx + cmprts->n_prts, d_id);
  cuda_mparticles_dump(cmprts);

  cuda_mparticles_reorder_and_offsets(cmprts);
  cuda_mparticles_dump(cmprts);
  
  cuda_mparticles_dealloc(cmprts);
#endif
  cuda_mparticles_destroy(cmprts);
}
