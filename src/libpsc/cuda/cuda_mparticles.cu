
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

cuda_mparticles::cuda_mparticles(const Grid_t& grid, const Int3& bs)
  : cuda_mparticles_base(grid, bs)
{
  xb_by_patch.resize(n_patches);
  for (int p = 0; p < n_patches; p++) {
    xb_by_patch[p] = grid.patches[p].xb;
  }

  cuda_mparticles_bnd::setup(this);
}

// ----------------------------------------------------------------------
// reserve_all

void cuda_mparticles::reserve_all(const uint *n_prts_by_patch)
{
  uint size = 0;
  for (int p = 0; p < n_patches; p++) {
    size += n_prts_by_patch[p];
  }

  if (size <= n_alloced) {
    return;
  }

  size *= 1.2;// FIXME hack
  n_alloced = std::max(size, 2 * n_alloced);

  cuda_mparticles_base::reserve_all();

  d_alt_xi4.resize(n_alloced);
  d_alt_pxi4.resize(n_alloced);
  d_bidx.resize(n_alloced);
  d_id.resize(n_alloced);

  cuda_mparticles_bnd::reserve_all(this);
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
      printf("cuda_mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g || bidx %d id %d\n",
	     n, xi4.x, xi4.y, xi4.z, cuda_float_as_int(xi4.w), pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	     bidx, id);
      assert(b == bidx);
    }
    off += off_e - off_b;
  }
}

// ----------------------------------------------------------------------
// cuda_mparticles_swap_alt

void
cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts)
{
  thrust::swap(cmprts->d_xi4, cmprts->d_alt_xi4);
  thrust::swap(cmprts->d_pxi4, cmprts->d_alt_pxi4);
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
    prm->b_mx[d]  = cuda_mprts->b_mx[d];
    prm->b_dxi[d] = cuda_mprts->b_dxi[d];
  }
}

static void
cuda_params2_free(struct cuda_params2 *prm)
{
}

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// get_block_idx

static int
get_block_idx(struct cuda_mparticles *cmprts, float4 xi4, int p)
{
  float *b_dxi = cmprts->b_dxi;
  int *b_mx = cmprts->b_mx;
  
  uint block_pos_y = (int) floorf(xi4.y * b_dxi[1]);
  uint block_pos_z = (int) floorf(xi4.z * b_dxi[2]);

  int bidx;
  if (block_pos_y >= b_mx[1] || block_pos_z >= b_mx[2]) {
    bidx = -1;
  } else {
    bidx = (p * b_mx[2] + block_pos_z) * b_mx[1] + block_pos_y;
  }

  return bidx;
}

// ----------------------------------------------------------------------
// cuda_mprts_find_block_indices_ids

__global__ static void
mprts_find_block_indices_ids(struct cuda_params2 prm, float4 *d_xi4, uint *d_off,
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

  struct cuda_params2 prm;
  cuda_params2_set(&prm, this);
    
  dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_find_block_indices_ids<<<dimGrid, dimBlock>>>(prm,
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
// cuda_mparticles_reorder_and_offsets

__global__ static void
mprts_reorder_and_offsets(int nr_prts, float4 *xi4, float4 *pxi4, float4 *alt_xi4, float4 *alt_pxi4,
			  uint *d_bidx, uint *d_ids, uint *d_off, int last_block)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i > nr_prts)
    return;

  int block, prev_block;
  if (i < nr_prts) {
    alt_xi4[i] = xi4[d_ids[i]];
    alt_pxi4[i] = pxi4[d_ids[i]];
    
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

void cuda_mparticles::reorder_and_offsets()
{
  if (n_patches == 0) {
    return;
  }

  dim3 dimGrid((n_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_reorder_and_offsets<<<dimGrid, dimBlock>>>(n_prts, d_xi4.data().get(), d_pxi4.data().get(),
						   d_alt_xi4.data().get(), d_alt_pxi4.data().get(),
						   d_bidx.data().get(), d_id.data().get(),
						   d_off.data().get(), n_blocks);
  cuda_sync_if_enabled();

  cuda_mparticles_swap_alt(this);
  need_reorder = false;
}

void
cuda_mparticles_reorder_and_offsets_slow(struct cuda_mparticles *cmprts)
{
  if (cmprts->n_patches == 0) {
    return;
  }

  thrust::host_vector<float4> h_xi4(cmprts->d_xi4.data(), cmprts->d_xi4.data() + cmprts->n_prts);
  thrust::host_vector<float4> h_pxi4(cmprts->d_pxi4.data(), cmprts->d_pxi4.data() + cmprts->n_prts);
  thrust::host_vector<float4> h_alt_xi4(cmprts->d_alt_xi4.data(), cmprts->d_alt_xi4.data() + cmprts->n_prts);
  thrust::host_vector<float4> h_alt_pxi4(cmprts->d_alt_pxi4.data(), cmprts->d_alt_pxi4.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_off(cmprts->d_off);
  thrust::host_vector<uint> h_bidx(cmprts->d_bidx.data(), cmprts->d_bidx.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_id(cmprts->d_id.data(), cmprts->d_id.data() + cmprts->n_prts);

  for (int i = 0; i <= cmprts->n_prts; i++) {
    //    uint bidx;
    uint block;
    if (i < cmprts->n_prts) {
      h_alt_xi4[i] = h_xi4[h_id[i]];
      h_alt_pxi4[i] = h_pxi4[h_id[i]];
      //bidx = get_block_idx(cmprts, h_alt_xi4[i], 0);
      block = h_bidx[i];
    } else {
      //bidx = cmprts->n_blocks;
      block = cmprts->n_blocks;
    }
    // if (i < 10) {
    //   printf("i %d bidx %d block %d xi4 %g %g\n", bidx, block, h_alt_xi4[i].y, h_alt_xi4[i].z);
    // }
    int prev_block = (i > 0) ? (int) h_bidx[i-1] : -1;
    for (int b = prev_block + 1; b <= block; b++) {
      h_off[b] = i;
    }
  }

  thrust::copy(h_alt_xi4.begin(), h_alt_xi4.end(), cmprts->d_alt_xi4.begin());
  thrust::copy(h_alt_pxi4.begin(), h_alt_pxi4.end(), cmprts->d_alt_pxi4.begin());
  thrust::copy(h_off.begin(), h_off.end(), cmprts->d_off.begin());
  
  cuda_mparticles_swap_alt(cmprts);
  cmprts->need_reorder = false;
}

// ----------------------------------------------------------------------
// check_in_patch_unordered_slow

void cuda_mparticles::check_in_patch_unordered_slow(uint *nr_prts_by_patch)
{
  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int n = 0; n < nr_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(this, d_xi4[off + n], p);
      assert(bidx >= 0 && bidx <= n_blocks);
    }
    off += nr_prts_by_patch[p];
  }

  assert(off == n_prts);
  printf("PASS: cuda_mparticles_check_in_patch_unordered_slow()\n");
}

// ----------------------------------------------------------------------
// check_bix_id_unordered_slow

void cuda_mparticles::check_bidx_id_unordered_slow(uint *n_prts_by_patch)
{
  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = get_block_idx(this, d_xi4[off + n], p);
      assert(bidx == d_bidx[off+n]);
      assert(off+n == d_id[off+n]);
    }
    off += n_prts_by_patch[p];
  }

  assert(off == n_prts);
  printf("PASS: cuda_mparticles_check_bidx_id_unordered_slow()\n");
}

// ----------------------------------------------------------------------
// check_ordered_slow

void cuda_mparticles::check_ordered_slow()
{
  uint off = 0;
  for (int b = 0; b < n_blocks; b++) {
    int p = b / n_blocks_per_patch;
    uint off_b = d_off[b], off_e = d_off[b+1];
    assert(off_e >= off_b);
    // printf("cuda_mparticles_check_ordered: block %d: %d -> %d (patch %d)\n", b, off_b, off_e, p);
    assert(d_off[b] == off);
    for (int n = d_off[b]; n < d_off[b+1]; n++) {
      float4 xi4;
      if (need_reorder) {
	xi4 = d_xi4[d_id[n]];
      } else {
	xi4 = d_xi4[n];
      }
      uint bidx = get_block_idx(this, xi4, p);
      //printf("cuda_mparticles_check_ordered: bidx %d\n", bidx);
      if (b != bidx) {
	printf("b %d bidx %d n %d p %d xi4 %g %g %g\n",
	       b, bidx, n, p, xi4.x, xi4.y, xi4.z);
	uint block_pos_y = (int) floorf(xi4.y * b_dxi[1]);
	uint block_pos_z = (int) floorf(xi4.z * b_dxi[2]);
	printf("block_pos %d %d %g %g\n", block_pos_y, block_pos_z, xi4.y * b_dxi[1],
	       xi4.z * b_dxi[2]);
      }
      assert(b == bidx);
    }
    off += off_e - off_b;
  }
  assert(off == n_prts);
  printf("cuda_mparticles_check_ordered: PASS\n");
}

// ----------------------------------------------------------------------
// check_ordered

void cuda_mparticles::check_ordered()
{
  thrust::host_vector<float4> h_xi4(d_xi4.data(), d_xi4.data() + n_prts);
  thrust::host_vector<uint> h_off(d_off);
  thrust::host_vector<uint> h_id(d_id.data(), d_id.data() + n_prts);

  //printf("cuda_mparticles_check_ordered: need_reorder %s\n", need_reorder ? "true" : "false");

  // for (int n = 0; n < 10; n++) {
  //   uint bidx = d_bidx[n];
  //   printf("n %d bidx %d xi4 %g %g\n", n, bidx, h_xi4[n].y, h_xi4[n].z);
  // }
  uint off = 0;
  for (int b = 0; b < n_blocks; b++) {
    int p = b / n_blocks_per_patch;
    uint off_b = h_off[b], off_e = h_off[b+1];
    assert(off_e >= off_b);
    //printf("cuda_mparticles_check_ordered: block %d: %d -> %d (patch %d)\n", b, off_b, off_e, p);
    assert(off_b == off);
    for (int n = h_off[b]; n < h_off[b+1]; n++) {
      float4 xi4;
      if (need_reorder) {
	xi4 = h_xi4[h_id[n]];
      } else {
	xi4 = h_xi4[n];
      }
      uint bidx = get_block_idx(this, xi4, p);
      //printf("cuda_mparticles_check_ordered: bidx %d\n", bidx);
      if (b != bidx) {
	printf("b %d bidx %d n %d p %d xi4 %g %g %g\n",
	       b, bidx, n, p, xi4.x, xi4.y, xi4.z);
	uint block_pos_y = (int) floorf(xi4.y * b_dxi[1]);
	uint block_pos_z = (int) floorf(xi4.z * b_dxi[2]);
	printf("block_pos %d %d %g %g\n", block_pos_y, block_pos_z, xi4.y * b_dxi[1],
	       xi4.z * b_dxi[2]);
      }
      assert(b == bidx);
    }
    off += off_e - off_b;
  }
  assert(off == n_prts);
  printf("cuda_mparticles_check_ordered: PASS\n");
}

// ----------------------------------------------------------------------
// cuda_mparticles_sort_initial

void
cuda_mparticles_sort_initial(struct cuda_mparticles *cmprts,
			     uint *n_prts_by_patch)
{
}

// ----------------------------------------------------------------------
// setup_internals

void cuda_mparticles::setup_internals()
{
  static int first_time = false;
  if (first_time) {
    uint n_prts_by_patch[n_patches];
    get_size_all(n_prts_by_patch);
    check_in_patch_unordered_slow(n_prts_by_patch);
  }

  find_block_indices_ids();
  if (first_time) {
    uint n_prts_by_patch[n_patches];
    get_size_all(n_prts_by_patch);
    check_bidx_id_unordered_slow(n_prts_by_patch);
  }

  thrust::stable_sort_by_key(d_bidx.data(), d_bidx.data() + n_prts, d_id.begin());
  reorder_and_offsets();

  if (first_time) {
    check_ordered();
    first_time = false;
  }
}

// ----------------------------------------------------------------------
// get_n_prts

uint cuda_mparticles::get_n_prts()
{
  return n_prts;
}

// ----------------------------------------------------------------------
// cuda_mparticles_reorder

__global__ static void
k_cuda_mparticles_reorder(int nr_prts, uint *d_ids,
		 float4 *xi4, float4 *pxi4,
		 float4 *alt_xi4, float4 *alt_pxi4)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i < nr_prts) {
    int j = d_ids[i];
    alt_xi4[i] = xi4[j];
    alt_pxi4[i] = pxi4[j];
  }
}

void
cuda_mparticles_reorder(struct cuda_mparticles *cmprts)
{
  if (!cmprts->need_reorder) {
    return;
  }
  
  dim3 dimGrid((cmprts->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  
  k_cuda_mparticles_reorder<<<dimGrid, THREADS_PER_BLOCK>>>
    (cmprts->n_prts, cmprts->d_id.data().get(),
     cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
     cmprts->d_alt_xi4.data().get(), cmprts->d_alt_pxi4.data().get());
  
  cuda_mparticles_swap_alt(cmprts);

  cmprts->need_reorder = false;
}

// ----------------------------------------------------------------------
// inject

void cuda_mparticles::inject(cuda_mparticles_prt *buf,
			     uint *buf_n_by_patch)
{
  if (need_reorder) {
    cuda_mparticles_reorder(this);
    need_reorder = false;
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
      cuda_mparticles_prt *prt = &buf[off + n];
      
      xi4->x  = prt->xi[0];
      xi4->y  = prt->xi[1];
      xi4->z  = prt->xi[2];
      xi4->w  = cuda_int_as_float(prt->kind);
      pxi4->x = prt->pxi[0];
      pxi4->y = prt->pxi[1];
      pxi4->z = prt->pxi[2];
      pxi4->w = prt->qni_wni;

      h_bidx[off + n] = get_block_idx(this, *xi4, p);
      h_id[off + n] = n_prts + off + n;
    }
    off += buf_n_by_patch[p];
  }
  assert(off == buf_n);

  uint n_prts_by_patch[n_patches];
  get_size_all(n_prts_by_patch);

  // check_in_patch_unordered_slow(n_prts_by_patch);

  find_block_indices_ids();
  // check_bidx_id_unordered_slow(n_prts_by_patch);

  assert(n_prts + buf_n <= n_alloced);
  thrust::copy(h_xi4.begin(), h_xi4.end(), d_xi4.data() + n_prts);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), d_pxi4.data() + n_prts);
  thrust::copy(h_bidx.begin(), h_bidx.end(), d_bidx.data() + n_prts);
  //thrust::copy(h_id.begin(), h_id.end(), d_id + n_prts);
  thrust::sequence(d_id.data(), d_id.data() + n_prts + buf_n);

  // for (int i = -5; i <= 5; i++) {
  //   //    float4 xi4 = d_xi4[cmprts->n_prts + i];
  //   uint bidx = d_bidx[cmprts->n_prts + i];
  //   uint id = d_id[cmprts->n_prts + i];
  //   printf("i %d bidx %d %d\n", i, bidx, id);
  // }

  // cmprts->check_ordered(cmprts);

  n_prts += buf_n;

  thrust::stable_sort_by_key(d_bidx.data(), d_bidx.data() + n_prts, d_id.begin());
  reorder_and_offsets();

  // cmprts->check_ordered(cmprts);
}

// ----------------------------------------------------------------------
// set_particles

void cuda_mparticles::set_particles(uint n_prts, uint off,
				    void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx),
				    void *ctx)
{
  float4 *xi4  = new float4[n_prts];
  float4 *pxi4 = new float4[n_prts];
  
  for (int n = 0; n < n_prts; n++) {
    struct cuda_mparticles_prt prt;
    get_particle(&prt, n, ctx);

    for (int d = 0; d < 3; d++) {
      int bi = fint(prt.xi[d] * b_dxi[d]);
      if (bi < 0 || bi >= b_mx[d]) {
	printf("XXX xi %g %g %g\n", prt.xi[0], prt.xi[1], prt.xi[2]);
	printf("XXX n %d d %d xi4[n] %g biy %d // %d\n",
	       n, d, prt.xi[d], bi, b_mx[d]);
	if (bi < 0) {
	  prt.xi[d] = 0.f;
	} else {
	  prt.xi[d] *= (1. - 1e-6);
	}
      }
      bi = floorf(prt.xi[d] * b_dxi[d]);
      assert(bi >= 0 && bi < b_mx[d]);
    }

    xi4[n].x  = prt.xi[0];
    xi4[n].y  = prt.xi[1];
    xi4[n].z  = prt.xi[2];
    xi4[n].w  = cuda_int_as_float(prt.kind);
    pxi4[n].x = prt.pxi[0];
    pxi4[n].y = prt.pxi[1];
    pxi4[n].z = prt.pxi[2];
    pxi4[n].w = prt.qni_wni;
  }

  to_device(xi4, pxi4, n_prts, off);
  
  delete[] xi4;
  delete[] pxi4;
}

// ----------------------------------------------------------------------
// get_particles

void cuda_mparticles::get_particles(uint n_prts, uint off,
				    void (*put_particle)(struct cuda_mparticles_prt *, int, void *),
				    void *ctx)
{
  float4 *xi4  = new float4[n_prts];
  float4 *pxi4 = new float4[n_prts];

  cuda_mparticles_reorder(this);
  from_device(xi4, pxi4, n_prts, off);
  
  for (int n = 0; n < n_prts; n++) {
    struct cuda_mparticles_prt prt;
    prt.xi[0]   = xi4[n].x;
    prt.xi[1]   = xi4[n].y;
    prt.xi[2]   = xi4[n].z;
    prt.kind    = cuda_float_as_int(xi4[n].w);
    prt.pxi[0]  = pxi4[n].x;
    prt.pxi[1]  = pxi4[n].y;
    prt.pxi[2]  = pxi4[n].z;
    prt.qni_wni = pxi4[n].w;

    put_particle(&prt, n, ctx);

#if 0
    for (int d = 0; d < 3; d++) {
      int bi = fint(prt.xi[d] * b_dxi[d]);
      if (bi < 0 || bi >= b_mx[d]) {
	MHERE;
	mprintf("XXX xi %.10g %.10g %.10g\n", prt.xi[0], prt.xi[1], prt.xi[2]);
	mprintf("XXX n %d d %d xi %.10g b_dxi %.10g bi %d // %d\n",
		n, d, prt.xi[d] * b_dxi[d], b_dxi[d], bi, b_mx[d]);
      }
    }
#endif
  }

  delete[] (xi4);
  delete[] (pxi4);
}

// ----------------------------------------------------------------------
// patch_get_b_dxi

const particle_cuda_real_t* cuda_mparticles::patch_get_b_dxi(int p)
{
  return b_dxi;
}

// ----------------------------------------------------------------------
// patch_get_b_mx

const int* cuda_mparticles::patch_get_b_mx(int p)
{
  return b_mx;
}

