
#include <cstdio>
#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "../libpsc/cuda/cuda_mparticles.h"

void
cuda_domain_info_set_test_1(struct cuda_domain_info *info)
{
  info->n_patches = 1;
  info->ldims[0] = 1; info->ldims[1] = 4; info->ldims[2] = 2;
  info->bs[0] = 1; info->bs[1] = 1; info->bs[2] = 1;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;
};

void
cuda_mparticles_add_particles_test_1(struct cuda_mparticles *cmprts,
				     unsigned int *n_prts_by_patch)
{
  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] = cmprts->ldims[0] * cmprts->ldims[1] * cmprts->ldims[2];
  }

  cuda_mparticles_reserve(cmprts, n_prts_by_patch);
  
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);

  int *ldims = cmprts->ldims;
  float *dx = cmprts->dx;
  
  int n = 0;
  for (int i = 0; i < ldims[0]; i++) {
    for (int j = 0; j < ldims[1]; j++) {
      for (int k = 0; k < ldims[2]; k++) {
	d_xi4[n] = (float4) { dx[0] * (i + .5f),
			      dx[1] * (j + .5f),
			      dx[1] * (k + .5f), 1. };
	d_pxi4[n] = (float4) { i, j, k, 2. };
	n++;
      }
    }
  }
}

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

int
main(void)
{
  struct cuda_mparticles *cmprts = cuda_mparticles_create();

  struct cuda_domain_info info;
  cuda_domain_info_set_test_1(&info);

  cuda_mparticles_set_domain_info(cmprts, &info);
  unsigned int n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_1(cmprts, n_prts_by_patch);
  printf("added particles\n");
  cuda_mparticles_dump(cmprts);
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
  cuda_mparticles_destroy(cmprts);
}
