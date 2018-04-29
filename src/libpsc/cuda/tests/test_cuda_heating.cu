
#include <cstdio>
#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include "../libpsc/cuda/cuda_mparticles.h"
#include "../libpsc/cuda/cuda_bits.h"

#include <mrc_profile.h>

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

// set up a domain [-40:40] x [-20:20], with 2 patches, cell size of 10

void
cuda_domain_info_set_test_2(struct cuda_domain_info *info)
{
  info->n_patches = 2;
  info->ldims[0] = 1; info->ldims[1] = 4; info->ldims[2] = 4;
  info->bs[0] = 1; info->bs[1] = 2; info->bs[2] = 2;
  info->dx[0] = 1.; info->dx[1] = 10.; info->dx[2] = 10.;

  info->xb_by_patch = new double_3[info->n_patches]; // FIXME, leaked
  info->xb_by_patch[0][0] = 0.;
  info->xb_by_patch[0][1] = -40.;
  info->xb_by_patch[0][2] = -20.;
  info->xb_by_patch[1][0] = 0.;
  info->xb_by_patch[1][1] = 0.;
  info->xb_by_patch[1][2] = -20.;
};

void
cuda_mparticles_add_particles_test_2(struct cuda_mparticles *cmprts,
				     uint *n_prts_by_patch)
{
  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] = 2 * cmprts->ldims[0] * cmprts->ldims[1] * cmprts->ldims[2];
  }

  cuda_mparticles_reserve_all(cmprts, n_prts_by_patch);
  
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);

  int *ldims = cmprts->ldims;
  float *dx = cmprts->dx;
  
  uint off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    uint n = 0;
    int ijk[3];
    for (ijk[0] = 0; ijk[0] < ldims[0]; ijk[0]++) {
      for (ijk[1] = 0; ijk[1] < ldims[1]; ijk[1]++) {
	for (ijk[2] = 0; ijk[2] < ldims[2]; ijk[2]++) {
	  double x[3];
	  for (int d = 0; d < 3; d++) {
	    x[d] = dx[d] * (ijk[d] + .5 );
	  }
	  for (int kind = 0; kind < 2; kind++) {
	    d_xi4[n + off] = (float4) { x[0], x[1], x[2], cuda_int_as_float(kind) };
	    d_pxi4[n + off] = (float4) { ijk[0], ijk[1], ijk[2], 99. };
	    n++;
	  }
	}
      }
    }
    assert(n == n_prts_by_patch[p]);
    off += n_prts_by_patch[p];
  }
}

static int
get_block_idx(struct cuda_mparticles *cmprts, int n, int p)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  int *b_mx = cmprts->b_mx;
  
  float4 xi4 = d_xi4[n];
  uint block_pos_y = cmprts->blockPosition(xi4.y, 1);
  uint block_pos_z = cmprts->blockPosition(xi4.z, 2);

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
					 uint *nr_prts_by_patch)
{
  uint off = 0;
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
					uint *n_prts_by_patch)
{
  thrust::device_ptr<uint> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<uint> d_id(cmprts->d_id);

  uint off = 0;
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

  struct cuda_domain_info info = {};
  cuda_domain_info_set_test_2(&info);

  cuda_mparticles_set_domain_info(cmprts, &info);
  uint n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_add_particles_test_2(cmprts, n_prts_by_patch);
  printf("added particles\n");
  cuda_mparticles_dump_by_patch(cmprts, n_prts_by_patch);
  cuda_mparticles_check_in_patch_unordered(cmprts, n_prts_by_patch);

  cuda_mparticles_find_block_indices_ids(cmprts);
  printf("find bidx, id\n");
  cuda_mparticles_dump_by_patch(cmprts, n_prts_by_patch);
  cuda_mparticles_check_bidx_id_unordered(cmprts, n_prts_by_patch);

  thrust::device_ptr<uint> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<uint> d_id(cmprts->d_id);
  thrust::stable_sort_by_key(d_bidx, d_bidx + cmprts->n_prts, d_id);
  printf("sort bidx, id\n");

  cuda_mparticles_reorder_and_offsets(cmprts);
  cuda_mparticles_dump(cmprts);

  // FIXME, up to here, it's been a copy of test_cuda_mparticles to set up
  // a valid particle structure

  HeatingSpotFoilParams foil_params;
  foil_params.zl = -10.f;
  foil_params.zh =  10.f;
  foil_params.xc = 0.f;
  foil_params.yc = 0.f;
  foil_params.rH = 30.f;
  foil_params.T  = .04f;
  foil_params.Mi = 100.f;

  printf("setup_foil\n");
  cuda_heating_foil foil(foil_params, 1, .1);

  printf("run_foil\n");
  cuda_heating_run_foil(cmprts);

  printf("after heating_run_foil\n");
  cuda_mparticles_dump(cmprts);

  cuda_mparticles_destroy(cmprts);
}
