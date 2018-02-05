
#include "cuda_mparticles.h"
#include "cuda_bits.h"

// ----------------------------------------------------------------------
// check_in_patch_unordered_slow

bool cuda_mparticles::check_in_patch_unordered_slow()
{
  uint n_prts_by_patch[n_patches];
  get_size_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(d_xi4[off + n], p);
      if (!(bidx >= 0 && bidx <= n_blocks)) return false;
    }
    off += n_prts_by_patch[p];
  }

  if (!(off == n_prts)) return false;
  // printf("PASS: cuda_mparticles_check_in_patch_unordered_slow()\n");
  return true;
}

// ----------------------------------------------------------------------
// check_bix_id_unordered_slow
//
// checks that block indices are correct,
// id is just enumerating particles

bool cuda_mparticles::check_bidx_id_unordered_slow()
{
  uint n_prts_by_patch[n_patches];
  get_size_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(d_xi4[off + n], p);
      if (!(bidx == d_bidx[off+n])) return false;
      if (!(off+n == d_id[off+n])) return false;
    }
    off += n_prts_by_patch[p];
  }

  if (!(off == n_prts)) return false;
  // printf("PASS: cuda_mparticles_check_bidx_id_unordered_slow()\n");
  return true;
}

// ----------------------------------------------------------------------
// check_ordered

bool cuda_mparticles::check_ordered()
{
  thrust::host_vector<float4> h_xi4(d_xi4.data(), d_xi4.data() + n_prts);
  thrust::host_vector<uint> h_off(d_off);
  thrust::host_vector<uint> h_id(d_id.data(), d_id.data() + n_prts);

  //printf("check_ordered: need_reorder %s\n", need_reorder ? "true" : "false");

  uint off = 0;
  for (int b = 0; b < n_blocks; b++) {
    int p = b / n_blocks_per_patch;
    uint off_b = h_off[b], off_e = h_off[b+1];
    if (!(off_e >= off_b)) return false;
    //printf("check_ordered: block %d: %d -> %d (patch %d)\n", b, off_b, off_e, p);
    if (!(off_b == off)) return false;
    for (int n = h_off[b]; n < h_off[b+1]; n++) {
      float4 xi4;
      if (need_reorder) {
	xi4 = h_xi4[h_id[n]];
      } else {
	xi4 = h_xi4[n];
      }
      uint bidx = get_block_idx(xi4, p);
      //printf("check_ordered: bidx %d\n", bidx);
      if (b != bidx) {
	printf("b %d bidx %d n %d p %d xi4 %g %g %g\n",
	       b, bidx, n, p, xi4.x, xi4.y, xi4.z);
	uint block_pos_y = blockPosition(xi4.y, 1);
	uint block_pos_z = blockPosition(xi4.z, 2);
	printf("block_pos %d %d %g %g\n", block_pos_y, block_pos_z,
	       xi4.y * pi_.b_dxi_[1], xi4.z * pi_.b_dxi_[2]);
      }
      if (!(b == bidx)) return false;
    }
    off += off_e - off_b;
  }
  if (!(off == n_prts)) return false;
  // printf("PASS: cuda_mparticles_check_ordered:\n");
  return true;
}

