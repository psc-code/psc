
#include <cstdio>
#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "cuda_mparticles.h"

void
cuda_domain_info_set_test_1(struct cuda_domain_info *info)
{
  info->nr_patches = 1;
  info->mx[0] = 1; info->mx[1] = 4; info->mx[2] = 2;
  info->bs[0] = 1; info->bs[1] = 1; info->bs[2] = 1;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;
  for (int d = 0; d < 3; d++) {
    assert(info->mx[d] % info->bs[d] == 0);
  }
};

void
cuda_mparticles_add_particles_test_1(struct cuda_mparticles *cuda_mprts,
				     unsigned int *nr_prts_by_patch)
{
  for (int p = 0; p < cuda_mprts->nr_patches; p++) {
    nr_prts_by_patch[p] = cuda_mprts->mx[0] * cuda_mprts->mx[1] * cuda_mprts->mx[2];
  }

  cuda_mparticles_alloc(cuda_mprts, nr_prts_by_patch);
  
  thrust::device_ptr<float4> d_xi4(cuda_mprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cuda_mprts->d_pxi4);

  int *mx = cuda_mprts->mx;
  float *dx = cuda_mprts->dx;
  
  int n = 0;
  for (int i = 0; i < mx[0]; i++) {
    for (int j = 0; j < mx[1]; j++) {
      for (int k = 0; k < mx[2]; k++) {
	d_xi4[n] = (float4) { dx[0] * (i + .5f),
			      dx[1] * (j + .5f),
			      dx[1] * (k + .5f), 1. };
	d_pxi4[n] = (float4) { i, j, k, 2. };
	n++;
      }
    }
  }
}

int
main(void)
{
  struct cuda_mparticles _cuda_mprts, *cuda_mprts = &_cuda_mprts;

  struct cuda_domain_info info;
  cuda_domain_info_set_test_1(&info);

  cuda_mparticles_set_domain_info(cuda_mprts, &info);
  unsigned int nr_prts_by_patch[cuda_mprts->nr_patches];
  cuda_mparticles_add_particles_test_1(cuda_mprts, nr_prts_by_patch);
  cuda_mparticles_dump(cuda_mprts);

  cuda_mparticles_find_block_indices_ids_total(cuda_mprts, nr_prts_by_patch);
  cuda_mparticles_dump(cuda_mprts);

  thrust::device_ptr<unsigned int> d_bidx(cuda_mprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cuda_mprts->d_id);
  thrust::stable_sort_by_key(d_bidx, d_bidx + cuda_mprts->nr_prts, d_id);
  cuda_mparticles_dump(cuda_mprts);

  cuda_mparticles_reorder_and_offsets(cuda_mprts);
  cuda_mparticles_dump(cuda_mprts);
  
  cuda_mparticles_free(cuda_mprts);
}
