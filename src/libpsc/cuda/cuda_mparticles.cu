
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

template<typename BS>
cuda_mparticles<BS>::cuda_mparticles(const Grid_t& grid)
: cuda_mparticles_base<BS>(grid)
{
  cuda_base_init();
  
  xb_by_patch.resize(this->n_patches);
  for (int p = 0; p < this->n_patches; p++) {
    xb_by_patch[p] = Real3(grid.patches[p].xb);
  }
}

// ----------------------------------------------------------------------
// reserve_all

template<typename BS>
void cuda_mparticles<BS>::reserve_all(const uint *n_prts_by_patch)
{
  uint size = 0;
  for (int p = 0; p < this->n_patches; p++) {
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

template<typename BS>
void cuda_mparticles<BS>::resize(uint n_prts)
{
  cuda_mparticles_base<BS>::reserve_all(n_prts);
  this->by_block_.d_idx.resize(n_prts);
  this->by_block_.d_id.resize(n_prts);
}

// ----------------------------------------------------------------------
// dump_by_patch

template<typename BS>
void cuda_mparticles<BS>::dump_by_patch(uint *n_prts_by_patch)
{
  printf("cuda_mparticles_dump_by_patch: n_prts = %d\n", this->n_prts);
  uint off = 0;
  for (int p = 0; p < this->n_patches; p++) {
    float *xb = &xb_by_patch[p][0];
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      float4 xi4 = this->d_xi4[n + off], pxi4 = this->d_pxi4[n + off];
      uint bidx = this->by_block_.d_idx[n + off], id = this->by_block_.d_id[n + off];
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

template<typename BS>
void cuda_mparticles<BS>::dump(const std::string& filename) const
{
  FILE* file = fopen(filename.c_str(), "w");
  assert(file);
  
  fprintf(file, "cuda_mparticles_dump: n_prts = %d\n", this->n_prts);
  uint off = 0;
  auto& d_off = this->by_block_.d_off;
  for (int b = 0; b < this->n_blocks; b++) {
    uint off_b = d_off[b], off_e = d_off[b+1];
    int p = b / this->n_blocks_per_patch;
    fprintf(file, "cuda_mparticles_dump: block %d: %d -> %d (patch %d)\n", b, off_b, off_e, p);
    assert(d_off[b] == off);
    for (int n = d_off[b]; n < d_off[b+1]; n++) {
      float4 xi4 = this->d_xi4[n], pxi4 = this->d_pxi4[n];
      uint bidx = this->by_block_.d_idx[n], id = this->by_block_.d_id[n];
      fprintf(file, "mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g || bidx %d id %d %s\n",
	      n, xi4.x, xi4.y, xi4.z, cuda_float_as_int(xi4.w), pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	      bidx, id, b == bidx ? "" : "BIDX MISMATCH!");
    }
    off += off_e - off_b;
  }
  fclose(file);
}

// ----------------------------------------------------------------------
// swap_alt

template<typename BS>
void cuda_mparticles<BS>::swap_alt()
{
  thrust::swap(this->d_xi4, d_alt_xi4);
  thrust::swap(this->d_pxi4, d_alt_pxi4);
}

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// k_find_block_indices_ids

template<typename BS>
__global__ static void
k_find_block_indices_ids(DParticleIndexer<BS> dpi, float4 *d_xi4, uint *d_off,
			 uint *d_bidx, uint *d_ids, int n_patches,
			 int n_blocks_per_patch)
{
  for (int p = 0; p < n_patches; p++) {
    uint off = d_off[p * n_blocks_per_patch];
    uint n_prts = d_off[(p + 1) * n_blocks_per_patch] - off;

    int n = threadIdx.x + blockDim.x * blockIdx.x;
    for (; n < n_prts; n += gridDim.x * blockDim.x) {
      float4 xi4 = d_xi4[n + off];
      d_bidx[n + off] = dpi.blockIndex(xi4, p);
      d_ids[n + off] = n + off;
    }
  }
}

// ----------------------------------------------------------------------
// k_find_cell_indices_ids

template<typename BS>
__global__ static void
k_find_cell_indices_ids(DParticleIndexer<BS> dpi, float4 *d_xi4, uint *d_off,
			uint *d_cidx, uint *d_ids, int n_patches,
			int n_blocks_per_patch)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  for (int p = 0; p < n_patches; p++) {
    uint off = d_off[p * n_blocks_per_patch];
    uint n_prts = d_off[(p + 1) * n_blocks_per_patch] - off;
    if (n < n_prts) {
      float4 xi4 = d_xi4[n + off];
      d_cidx[n + off] = dpi.validCellIndex(xi4, p);
      d_ids[n + off] = n + off;
    }
  }
}

// ----------------------------------------------------------------------
// find_block_indices_ids

template<typename BS>
void cuda_mparticles<BS>::find_block_indices_ids(thrust::device_vector<uint>& d_idx,
						 thrust::device_vector<uint>& d_id)
{
  if (this->n_patches == 0) {
    return;
  }

  // OPT: if we didn't need max_n_prts, we wouldn't have to get the
  // sizes / offsets at all, and it seems likely we could do a better
  // job here in general
  uint n_prts_by_patch[this->n_patches];
  this->get_size_all(n_prts_by_patch);
  
  int max_n_prts = 0;
  for (int p = 0; p < this->n_patches; p++) {
    if (n_prts_by_patch[p] > max_n_prts) {
      max_n_prts = n_prts_by_patch[p];
    }
  }

  if (max_n_prts == 0) {
    return;
  }

  int n_blocks = (max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  if (n_blocks > 32768) n_blocks = 32768;
  dim3 dimGrid(n_blocks);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_find_block_indices_ids<BS><<<dimGrid, dimBlock>>>(*this,
						      this->d_xi4.data().get(),
						      this->by_block_.d_off.data().get(),
						      d_idx.data().get(),
						      d_id.data().get(),
						      this->n_patches,
						      this->n_blocks_per_patch);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// find_cell_indices_ids

template<typename BS>
void cuda_mparticles<BS>::find_cell_indices_ids(thrust::device_vector<uint>& d_cidx,
						thrust::device_vector<uint>& d_id)
{
  if (this->n_patches == 0) {
    return;
  }

  // OPT: if we didn't need max_n_prts, we wouldn't have to get the
  // sizes / offsets at all, and it seems likely we could do a better
  // job here in general
  uint n_prts_by_patch[this->n_patches];
  this->get_size_all(n_prts_by_patch);
  
  int max_n_prts = 0;
  for (int p = 0; p < this->n_patches; p++) {
    if (n_prts_by_patch[p] > max_n_prts) {
      max_n_prts = n_prts_by_patch[p];
    }
  }

  if (max_n_prts == 0) {
    return;
  }
  
  dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_find_cell_indices_ids<BS><<<dimGrid, dimBlock>>>(*this,
						     this->d_xi4.data().get(),
						     this->by_block_.d_off.data().get(),
						     d_cidx.data().get(),
						     d_id.data().get(),
						     this->n_patches,
						     this->n_blocks_per_patch);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// k_reorder_and_offsets

__global__ static void
k_reorder_and_offsets(int nr_prts, float4 *xi4, float4 *pxi4, float4 *alt_xi4, float4 *alt_pxi4,
		      const uint *d_bidx, const uint *d_ids, uint *d_off, int last_block)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;

  for (; i <= nr_prts; i += blockDim.x * gridDim.x) {
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
}

// ----------------------------------------------------------------------
// reorder_and_offsets

template<typename BS>
void cuda_mparticles<BS>::reorder_and_offsets(const thrust::device_vector<uint>& d_idx,
					      const thrust::device_vector<uint>& d_id,
					      thrust::device_vector<uint>& d_off)
{
  if (this->n_patches == 0) {
    return;
  }

  swap_alt();
  resize(this->n_prts);

  int n_blocks = (this->n_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  if (n_blocks > 32768) n_blocks = 32768;
  dim3 dimGrid(n_blocks);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_reorder_and_offsets<<<dimGrid, dimBlock>>>(this->n_prts, this->d_xi4.data().get(), this->d_pxi4.data().get(),
					       d_alt_xi4.data().get(), d_alt_pxi4.data().get(),
					       d_idx.data().get(),
					       d_id.data().get(),
					       d_off.data().get(), this->n_blocks);
  cuda_sync_if_enabled();

  need_reorder = false;
}

// ----------------------------------------------------------------------
// k_reorder

__global__ static void
k_reorder(int n_prts, const uint *d_ids, float4 *xi4, float4 *pxi4,
	  const float4 *alt_xi4, const float4 *alt_pxi4)
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

template<typename BS>
void cuda_mparticles<BS>::reorder()
{
  if (!need_reorder) {
    return;
  }

  reorder(this->by_block_.d_id);
  need_reorder = false;
}

// ----------------------------------------------------------------------
// reorder

template<typename BS>
void cuda_mparticles<BS>::reorder(const thrust::device_vector<uint>& d_id)
{
  swap_alt();
  resize(this->n_prts);

  dim3 dimGrid((this->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  
  k_reorder<<<dimGrid, THREADS_PER_BLOCK>>>
    (this->n_prts, d_id.data().get(), this->d_xi4.data().get(), this->d_pxi4.data().get(),
     d_alt_xi4.data().get(), d_alt_pxi4.data().get());
}

// ----------------------------------------------------------------------
// setup_internals

template<typename BS>
void cuda_mparticles<BS>::setup_internals()
{
  // pre-condition: particles sorted by patch, d_off being used to
  // describe patch boundaries

  // assert(check_in_patch_unordered_slow());

  this->by_block_.find_indices_ids(*this);

  // assert(check_bidx_id_unordered_slow());

  this->by_block_.stable_sort();

  this->by_block_.reorder_and_offsets(*this);

  // post-condition:
  // - particles now sorted by block
  // - d_off describes block boundaries
  // - UNUSED: d_bidx has each particle's block index

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// get_n_prts

template<typename BS>
uint cuda_mparticles<BS>::get_n_prts()
{
  return this->n_prts;
}

// ----------------------------------------------------------------------
// inject_buf

template<typename BS>
void cuda_mparticles<BS>::inject_buf(const cuda_mparticles_prt *buf,
				     const uint *buf_n_by_patch)
{
  if (need_reorder) {
    reorder();
  }
  
  uint buf_n = 0;
  for (int p = 0; p < this->n_patches; p++) {
    buf_n += buf_n_by_patch[p];
    //    printf("p %d buf_n_by_patch %d\n", p, buf_n_by_patch[p]);
  }
  //  printf("buf_n %d\n", buf_n);

  thrust::host_vector<float4> h_xi4(buf_n);
  thrust::host_vector<float4> h_pxi4(buf_n);
  thrust::host_vector<uint> h_bidx(buf_n);
  thrust::host_vector<uint> h_id(buf_n);

  uint off = 0;
  for (int p = 0; p < this->n_patches; p++) {
    for (int n = 0; n < buf_n_by_patch[p]; n++) {
      float4 *xi4 = &h_xi4[off + n];
      float4 *pxi4 = &h_pxi4[off + n];
      const cuda_mparticles_prt *prt = &buf[off + n];
      
      xi4->x  = prt->x[0];
      xi4->y  = prt->x[1];
      xi4->z  = prt->x[2];
      xi4->w  = cuda_int_as_float(prt->kind);
      pxi4->x = prt->p[0];
      pxi4->y = prt->p[1];
      pxi4->z = prt->p[2];
      pxi4->w = prt->w * this->grid_.kinds[prt->kind].q;

      auto bidx = this->blockIndex(*xi4, p);
      assert(bidx >= 0 && bidx < this->n_blocks);
      h_bidx[off + n] = bidx;;
      h_id[off + n] = this->n_prts + off + n;
    }
    off += buf_n_by_patch[p];
  }
  assert(off == buf_n);

  // assert(check_in_patch_unordered_slow());

  this->by_block_.find_indices_ids(*this);
  // assert(check_bidx_id_unordered_slow());

  resize(this->n_prts + buf_n);

  thrust::copy(h_xi4.begin(), h_xi4.end(), this->d_xi4.begin() + this->n_prts);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), this->d_pxi4.begin() + this->n_prts);
  thrust::copy(h_bidx.begin(), h_bidx.end(), this->by_block_.d_idx.begin() + this->n_prts);
  //thrust::copy(h_id.begin(), h_id.end(), d_id + n_prts);
  // FIXME, looks like ids up until n_prts have already been set above
  thrust::sequence(this->by_block_.d_id.data(), this->by_block_.d_id.data() + this->n_prts + buf_n);

  // for (int i = -5; i <= 5; i++) {
  //   //    float4 xi4 = d_xi4[cmprts->n_prts + i];
  //   uint bidx = d_bidx[cmprts->n_prts + i];
  //   uint id = d_id[cmprts->n_prts + i];
  //   printf("i %d bidx %d %d\n", i, bidx, id);
  // }

  // assert(check_ordered());

  this->n_prts += buf_n;

  this->by_block_.stable_sort();

  this->by_block_.reorder_and_offsets(*this);

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// inject_buf

template<typename BS>
void cuda_mparticles<BS>::inject_buf(const particle_inject *buf,
				     const uint *buf_n_by_patch)
{
  using Double3 = Vec3<double>;
  
  if (need_reorder) {
    reorder();
  }
  
  uint buf_n = 0;
  for (int p = 0; p < this->n_patches; p++) {
    buf_n += buf_n_by_patch[p];
    //    printf("p %d buf_n_by_patch %d\n", p, buf_n_by_patch[p]);
  }
  //  printf("buf_n %d\n", buf_n);

  thrust::host_vector<float4> h_xi4(buf_n);
  thrust::host_vector<float4> h_pxi4(buf_n);
  thrust::host_vector<uint> h_bidx(buf_n);
  thrust::host_vector<uint> h_id(buf_n);

  uint off = 0;
  for (int p = 0; p < this->n_patches; p++) {
    auto& patch = this->grid_.patches[p];
    for (int n = 0; n < buf_n_by_patch[p]; n++) {
      float4 *xi4 = &h_xi4[off + n];
      float4 *pxi4 = &h_pxi4[off + n];
      auto new_prt = buf[off + n];
      auto x = Double3{new_prt.x} - patch.xb;
      auto prt = cuda_mparticles_prt{Real3{x}, Real3{Double3{new_prt.u}},
				     real_t(new_prt.w), new_prt.kind};
  
      xi4->x  = prt.x[0];
      xi4->y  = prt.x[1];
      xi4->z  = prt.x[2];
      xi4->w  = cuda_int_as_float(prt.kind);
      pxi4->x = prt.p[0];
      pxi4->y = prt.p[1];
      pxi4->z = prt.p[2];
      pxi4->w = prt.w * this->grid_.kinds[prt.kind].q;

      auto bidx = this->blockIndex(*xi4, p);
      assert(bidx >= 0 && bidx < this->n_blocks);
      h_bidx[off + n] = bidx;;
      h_id[off + n] = this->n_prts + off + n;
    }
    off += buf_n_by_patch[p];
  }
  assert(off == buf_n);

  // assert(check_in_patch_unordered_slow());

  this->by_block_.find_indices_ids(*this);
  // assert(check_bidx_id_unordered_slow());

  resize(this->n_prts + buf_n);

  thrust::copy(h_xi4.begin(), h_xi4.end(), this->d_xi4.begin() + this->n_prts);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), this->d_pxi4.begin() + this->n_prts);
  thrust::copy(h_bidx.begin(), h_bidx.end(), this->by_block_.d_idx.begin() + this->n_prts);
  //thrust::copy(h_id.begin(), h_id.end(), d_id + n_prts);
  // FIXME, looks like ids up until n_prts have already been set above
  thrust::sequence(this->by_block_.d_id.data(), this->by_block_.d_id.data() + this->n_prts + buf_n);

  // for (int i = -5; i <= 5; i++) {
  //   //    float4 xi4 = d_xi4[cmprts->n_prts + i];
  //   uint bidx = d_bidx[cmprts->n_prts + i];
  //   uint id = d_id[cmprts->n_prts + i];
  //   printf("i %d bidx %d %d\n", i, bidx, id);
  // }

  // assert(check_ordered());

  this->n_prts += buf_n;

  this->by_block_.stable_sort();

  this->by_block_.reorder_and_offsets(*this);

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// get_particles

template<typename BS>
std::vector<cuda_mparticles_prt> cuda_mparticles<BS>::get_particles(int beg, int end)
{
  int n_prts = end - beg;
  std::vector<cuda_mparticles_prt> prts;
  prts.reserve(n_prts);

  reorder(); // FIXME? by means of this, this function disturbs the state...

  thrust::host_vector<float4> xi4(&this->d_xi4[beg], &this->d_xi4[end]);
  thrust::host_vector<float4> pxi4(&this->d_pxi4[beg], &this->d_pxi4[end]);

  for (int n = 0; n < n_prts; n++) {
    int kind = cuda_float_as_int(xi4[n].w);
    prts.emplace_back(Real3{xi4[n].x, xi4[n].y, xi4[n].z},
	              Real3{pxi4[n].x, pxi4[n].y, pxi4[n].z},
	              pxi4[n].w / float(this->grid_.kinds[kind].q), kind);

#if 0
    uint b = blockIndex(xi4[n], p);
    assert(b < n_blocks);
#endif
  }

  return prts;
}

// ----------------------------------------------------------------------
// get_particles

template<typename BS>
std::vector<cuda_mparticles_prt> cuda_mparticles<BS>::get_particles(int p)
{
  // FIXME, doing the copy here all the time would be nice to avoid
  // making sure we actually have a valid d_off would't hurt, either
  thrust::host_vector<uint> h_off(this->by_block_.d_off);

  uint beg = h_off[p * this->n_blocks_per_patch];
  uint end = h_off[(p+1) * this->n_blocks_per_patch];

  return get_particles(beg, end);
}

// ----------------------------------------------------------------------
// start

template<typename BS>
uint cuda_mparticles<BS>::start(int p)
{
  return this->by_block_.d_off[p * this->n_blocks_per_patch];
}

#include "cuda_mparticles_gold.cu"
#include "cuda_mparticles_checks.cu"

template struct cuda_mparticles<BS144>;
template struct cuda_mparticles<BS444>;
