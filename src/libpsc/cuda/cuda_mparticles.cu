
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include "psc_bits.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// ctor

cuda_mparticles::cuda_mparticles(const Grid_t& grid)
  : cuda_mparticles_base(grid)
{
  xb_by_patch.resize(n_patches);
  for (int p = 0; p < n_patches; p++) {
    xb_by_patch[p] = Real3(grid.patches[p].xb);
  }
}

// ----------------------------------------------------------------------
// reserve_all

void cuda_mparticles::reserve_all(const uint *n_prts_by_patch)
{
  uint size = 0;
  for (int p = 0; p < n_patches; p++) {
    size += n_prts_by_patch[p];
  }

  // FIXME, arguably this should just reserve here, not actually resize
  resize(size);
}

// ----------------------------------------------------------------------
// resize
//
// the goal here is to have d_xi4, d_pxi4, d_bidx and d_id always
// have the same size.

void cuda_mparticles::resize(uint n_prts)
{
  cuda_mparticles_base::reserve_all(n_prts);
  d_bidx.resize(n_prts);
  d_id.resize(n_prts);
}

// ----------------------------------------------------------------------
// dump_by_patch

void cuda_mparticles::dump_by_patch(uint *n_prts_by_patch)
{
  printf("cuda_mparticles_dump_by_patch: n_prts = %d\n", n_prts);
  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    float *xb = &xb_by_patch[p][0];
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      float4 xi4 = d_xi4[n + off], pxi4 = d_pxi4[n + off];
      uint bidx = d_bidx[n + off], id = d_id[n + off];
      printf("cuda_mparticles_dump_by_patch: [%d/%d] %g %g %g // %d // %g %g %g // %g b_idx %d id %d\n",
	     p, n, xi4.x + xb[0], xi4.y + xb[1], xi4.z + xb[2],
	     cuda_float_as_int(xi4.w),
	     pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	     bidx, id);
    }
    off += n_prts_by_patch[p];
  }
}

// ----------------------------------------------------------------------
// dump

void cuda_mparticles::dump()
{
  printf("cuda_mparticles_dump: n_prts = %d\n", n_prts);
  uint off = 0;
  for (int b = 0; b < n_blocks; b++) {
    uint off_b = d_off[b], off_e = d_off[b+1];
    int p = b / n_blocks_per_patch;
    printf("cuda_mparticles_dump: block %d: %d -> %d (patch %d)\n", b, off_b, off_e, p);
    assert(d_off[b] == off);
    for (int n = d_off[b]; n < d_off[b+1]; n++) {
      float4 xi4 = d_xi4[n], pxi4 = d_pxi4[n];
      uint bidx = d_bidx[n], id = d_id[n];
      printf("cuda_mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g || bidx %d id %d %s\n",
	     n, xi4.x, xi4.y, xi4.z, cuda_float_as_int(xi4.w), pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	     bidx, id, b == bidx ? "" : "BIDX MISMATCH!");
    }
    off += off_e - off_b;
  }
}

// ----------------------------------------------------------------------
// swap_alt

void cuda_mparticles::swap_alt()
{
  thrust::swap(d_xi4, d_alt_xi4);
  thrust::swap(d_pxi4, d_alt_pxi4);
}

// ----------------------------------------------------------------------
// cuda_params2

struct cuda_params2 {
  uint b_mx[3];
  float b_dxi[3];
};

static void
cuda_params2_set(struct cuda_params2 *prm, const struct cuda_mparticles *cuda_mprts)
{
  for (int d = 0; d < 3; d++) {
    prm->b_mx[d]  = cuda_mprts->pi_.b_mx_[d];
    prm->b_dxi[d] = cuda_mprts->pi_.b_dxi_[d];
  }
}

static void
cuda_params2_free(struct cuda_params2 *prm)
{
}

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// k_find_block_indices_ids

__global__ static void
k_find_block_indices_ids(struct cuda_params2 prm, float4 *d_xi4, uint *d_off,
			 uint *d_bidx, uint *d_ids, int n_patches,
			 int n_blocks_per_patch)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  for (int p = 0; p < n_patches; p++) {
    uint off = d_off[p * n_blocks_per_patch];
    uint n_prts = d_off[(p + 1) * n_blocks_per_patch] - off;
    if (n < n_prts) {
      float4 xi4 = d_xi4[n + off];
      uint block_pos_y = __float2int_rd(xi4.y * prm.b_dxi[1]);
      uint block_pos_z = __float2int_rd(xi4.z * prm.b_dxi[2]);
      
      int block_idx;
      if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
	block_idx = -1; // not supposed to happen here!
      } else {
	block_idx = block_pos_z * prm.b_mx[1] + block_pos_y + p * nr_blocks;
      }
      d_bidx[n + off] = block_idx;
      d_ids[n + off] = n + off;
    }
  }
}

// ----------------------------------------------------------------------
// find_block_indices_ids

void cuda_mparticles::find_block_indices_ids()
{
  if (n_patches == 0) {
    return;
  }

  // OPT: if we didn't need max_n_prts, we wouldn't have to get the
  // sizes / offsets at all, and it seems likely we could do a better
  // job here in general
  uint n_prts_by_patch[n_patches];
  get_size_all(n_prts_by_patch);
  
  int max_n_prts = 0;
  for (int p = 0; p < n_patches; p++) {
    if (n_prts_by_patch[p] > max_n_prts) {
      max_n_prts = n_prts_by_patch[p];
    }
  }

  if (max_n_prts == 0) {
    return;
  }
  
  struct cuda_params2 prm;
  cuda_params2_set(&prm, this);
    
  dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_find_block_indices_ids<<<dimGrid, dimBlock>>>(prm,
						  d_xi4.data().get(),
						  d_off.data().get(),
						  d_bidx.data().get(),
						  d_id.data().get(),
						  n_patches,
						  n_blocks_per_patch);
  cuda_sync_if_enabled();
  cuda_params2_free(&prm);
}

// ----------------------------------------------------------------------
// stable_sort_by_key

void cuda_mparticles::stable_sort_by_key()
{
  thrust::stable_sort_by_key(d_bidx.data(), d_bidx.data() + n_prts, d_id.begin());
}

// ----------------------------------------------------------------------
// k_reorder_and_offsets

__global__ static void
k_reorder_and_offsets(int nr_prts, float4 *xi4, float4 *pxi4, float4 *alt_xi4, float4 *alt_pxi4,
		      uint *d_bidx, uint *d_ids, uint *d_off, int last_block)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i > nr_prts)
    return;

  int block, prev_block;
  if (i < nr_prts) {
    xi4[i] = alt_xi4[d_ids[i]];
    pxi4[i] = alt_pxi4[d_ids[i]];
    
    block = d_bidx[i];
  } else { // needed if there is no particle in the last block
    block = last_block;
  }

  // OPT: d_bidx[i-1] could use shmem
  // create offsets per block into particle array
  prev_block = -1;
  if (i > 0) {
    prev_block = d_bidx[i-1];
  }
  for (int b = prev_block + 1; b <= block; b++) {
    d_off[b] = i;
  }
}

// ----------------------------------------------------------------------
// reorder_and_offsets

void cuda_mparticles::reorder_and_offsets()
{
  if (n_patches == 0) {
    return;
  }

  swap_alt();
  resize(n_prts);

  dim3 dimGrid((n_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_reorder_and_offsets<<<dimGrid, dimBlock>>>(n_prts, d_xi4.data().get(), d_pxi4.data().get(),
					       d_alt_xi4.data().get(), d_alt_pxi4.data().get(),
					       d_bidx.data().get(), d_id.data().get(),
					       d_off.data().get(), n_blocks);
  cuda_sync_if_enabled();

  need_reorder = false;
}

// ----------------------------------------------------------------------
// k_reorder

__global__ static void
k_reorder(int n_prts, uint *d_ids, float4 *xi4, float4 *pxi4,
	  float4 *alt_xi4, float4 *alt_pxi4)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i < n_prts) {
    int j = d_ids[i];
    xi4[i] = alt_xi4[j];
    pxi4[i] = alt_pxi4[j];
  }
}

// ----------------------------------------------------------------------
// reorder

void cuda_mparticles::reorder()
{
  if (!need_reorder) {
    return;
  }
  
  swap_alt();
  resize(n_prts);

  dim3 dimGrid((n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  
  k_reorder<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_prts, d_id.data().get(), d_xi4.data().get(), d_pxi4.data().get(),
     d_alt_xi4.data().get(), d_alt_pxi4.data().get());
  
  need_reorder = false;
}

// ----------------------------------------------------------------------
// setup_internals

void cuda_mparticles::setup_internals()
{
  // pre-condition: particles sorted by patch, d_off being used to
  // describe patch boundaries

  // assert(check_in_patch_unordered_slow());

  find_block_indices_ids();

  // assert(check_bidx_id_unordered_slow());

  stable_sort_by_key();

  reorder_and_offsets();

  // post-condition:
  // - particles now sorted by block
  // - d_off describes block boundaries
  // - UNUSED: d_bidx has each particle's block index

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// get_n_prts

uint cuda_mparticles::get_n_prts()
{
  return n_prts;
}

// ----------------------------------------------------------------------
// inject

void cuda_mparticles::inject(const cuda_mparticles_prt *buf,
			     uint *buf_n_by_patch)
{
  if (need_reorder) {
    reorder();
  }
  
  uint buf_n = 0;
  for (int p = 0; p < n_patches; p++) {
    buf_n += buf_n_by_patch[p];
    //    printf("p %d buf_n_by_patch %d\n", p, buf_n_by_patch[p]);
  }
  //  printf("buf_n %d\n", buf_n);

  thrust::host_vector<float4> h_xi4(buf_n);
  thrust::host_vector<float4> h_pxi4(buf_n);
  thrust::host_vector<uint> h_bidx(buf_n);
  thrust::host_vector<uint> h_id(buf_n);

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int n = 0; n < buf_n_by_patch[p]; n++) {
      float4 *xi4 = &h_xi4[off + n];
      float4 *pxi4 = &h_pxi4[off + n];
      const cuda_mparticles_prt *prt = &buf[off + n];
      
      xi4->x  = prt->xi[0];
      xi4->y  = prt->xi[1];
      xi4->z  = prt->xi[2];
      xi4->w  = cuda_int_as_float(prt->kind);
      pxi4->x = prt->pxi[0];
      pxi4->y = prt->pxi[1];
      pxi4->z = prt->pxi[2];
      pxi4->w = prt->qni_wni;

      h_bidx[off + n] = blockIndex(*xi4, p);
      h_id[off + n] = n_prts + off + n;
    }
    off += buf_n_by_patch[p];
  }
  assert(off == buf_n);

  // assert(check_in_patch_unordered_slow());

  find_block_indices_ids();
  // assert(check_bidx_id_unordered_slow());

  resize(n_prts + buf_n);

  thrust::copy(h_xi4.begin(), h_xi4.end(), d_xi4.begin() + n_prts);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), d_pxi4.begin() + n_prts);
  thrust::copy(h_bidx.begin(), h_bidx.end(), d_bidx.begin() + n_prts);
  //thrust::copy(h_id.begin(), h_id.end(), d_id + n_prts);
  // FIXME, looks like ids up until n_prts have already been set above
  thrust::sequence(d_id.data(), d_id.data() + n_prts + buf_n);

  // for (int i = -5; i <= 5; i++) {
  //   //    float4 xi4 = d_xi4[cmprts->n_prts + i];
  //   uint bidx = d_bidx[cmprts->n_prts + i];
  //   uint id = d_id[cmprts->n_prts + i];
  //   printf("i %d bidx %d %d\n", i, bidx, id);
  // }

  // assert(check_ordered());

  n_prts += buf_n;

  stable_sort_by_key();

  reorder_and_offsets();

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// patch_get_b_mx

const int* cuda_mparticles::patch_get_b_mx(int p)
{
  return pi_.b_mx_;
}

// ----------------------------------------------------------------------
// cast to DMParticles

cuda_mparticles::operator DMParticles()
{
  return DMParticles(d_xi4.data().get(), d_pxi4.data().get(),
		     d_alt_xi4.data().get(), d_alt_pxi4.data().get(),
		     d_off.data().get(), d_bidx.data().get(),
		     d_id.data().get(), n_blocks);
}


#include "cuda_mparticles_gold.cu"
#include "cuda_mparticles_checks.cu"