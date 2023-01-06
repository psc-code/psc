
#include "cuda_mparticles.hxx"
#include "cuda_bits.h"

#include "psc_bits.h"
#include "bs.hxx"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h>

#include "cuda_base.hxx"

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// ctor

template <typename BS>
cuda_mparticles<BS>::cuda_mparticles(const Grid_t& grid) : base_type(grid)
{
  cuda_base_init();

  xb_by_patch.resize(this->n_patches());
  for (int p = 0; p < this->n_patches(); p++) {
    xb_by_patch[p] = Real3(grid.patches[p].xb);
  }
}

template <typename BS>
cuda_mparticles<BS>::~cuda_mparticles()
{
  mem_particles -= allocated_bytes(this->by_block_.d_idx);
  mem_particles -= allocated_bytes(this->by_block_.d_id);
}

// ----------------------------------------------------------------------
// resize
//
// the goal here is to have d_xi4, d_pxi4, d_bidx and d_id always
// have the same size.

template <typename BS>
void cuda_mparticles<BS>::resize(uint n_prts)
{
  base_type::resize(n_prts);
  mem_particles -= allocated_bytes(this->by_block_.d_idx);
  this->by_block_.d_idx.resize(n_prts);
  mem_particles += allocated_bytes(this->by_block_.d_idx);
  mem_particles -= allocated_bytes(this->by_block_.d_id);
  this->by_block_.d_id.resize(n_prts);
  mem_particles += allocated_bytes(this->by_block_.d_id);
}

// ----------------------------------------------------------------------
// dump_by_patch

template <typename BS>
void cuda_mparticles<BS>::dump_by_patch(uint* n_prts_by_patch)
{
  printf("cuda_mparticles_dump_by_patch: n_prts = %d\n", this->n_prts);
  uint off = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    float* xb = &xb_by_patch[p][0];
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      auto prt = this->storage[n + off];
      uint bidx = this->by_block_.d_idx[n + off],
           id = this->by_block_.d_id[n + off];
      printf("cuda_mparticles_dump_by_patch: [%d/%d] %g %g %g // %d // %g %g "
             "%g // %g b_idx %d id %d\n",
             p, n, prt.x[0] + xb[0], prt.x[1] + xb[1], prt.x[2] + xb[2],
             prt.kind, prt.u[0], prt.u[1], prt.u[2], prt.qni_wni, bidx, id);
    }
    off += n_prts_by_patch[p];
  }
}

// ----------------------------------------------------------------------
// dump

template <typename BS>
void cuda_mparticles<BS>::dump(const std::string& filename) const
{
  FILE* file = fopen(filename.c_str(), "w");
  assert(file);

  fprintf(file, "cuda_mparticles_dump: n_prts = %d\n", this->n_prts);
  uint off = 0;
  auto& d_off = this->by_block_.d_off;
  for (int b = 0; b < this->n_blocks; b++) {
    uint off_b = d_off[b], off_e = d_off[b + 1];
    int p = b / this->n_blocks_per_patch;
    fprintf(file, "cuda_mparticles_dump: block %d: %d -> %d (patch %d)\n", b,
            off_b, off_e, p);
    assert(d_off[b] == off);
    for (int n = d_off[b]; n < d_off[b + 1]; n++) {
      auto prt = this->storage[n + off];
      uint bidx = this->by_block_.d_idx[n], id = this->by_block_.d_id[n];
      fprintf(file,
              "mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g || bidx "
              "%d id %d %s\n",
              n, prt.x[0], prt.x[1], prt.x[2], prt.kind, prt.u[0], prt.u[1],
              prt.u[2], prt.qni_wni, bidx, id,
              b == bidx ? "" : "BIDX MISMATCH!");
    }
    off += off_e - off_b;
  }
  fclose(file);
}

// ----------------------------------------------------------------------
// swap_alt

template <typename BS>
void cuda_mparticles<BS>::swap_alt()
{
  using std::swap;
  swap(this->storage, this->alt_storage);
}

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// k_reorder_and_offsets

template <typename BS>
__global__ static void k_reorder_and_offsets(DMparticlesCuda<BS> dmprts,
                                             int nr_prts, const uint* d_bidx,
                                             const uint* d_ids, uint* d_off,
                                             int last_block)
{
  int i = threadIdx.x + blockDim.x * blockIdx.x;

  for (; i <= nr_prts; i += blockDim.x * gridDim.x) {
    int block, prev_block;
    if (i < nr_prts) {
      dmprts.storage.store(dmprts.alt_storage[d_ids[i]],
                           i); // storage[i] = alt_storage[d_ids[i]];

      block = d_bidx[i];
    } else { // needed if there is no particle in the last block
      block = last_block;
    }

    // OPT: d_bidx[i-1] could use shmem
    // create offsets per block into particle array
    prev_block = -1;
    if (i > 0) {
      prev_block = d_bidx[i - 1];
    }
    for (int b = prev_block + 1; b <= block; b++) {
      d_off[b] = i;
    }
  }
}

// ----------------------------------------------------------------------
// reorder_and_offsets

template <typename BS>
void cuda_mparticles<BS>::reorder_and_offsets(
  const psc::device_vector<uint>& d_idx, const psc::device_vector<uint>& d_id,
  psc::device_vector<uint>& d_off)
{
  if (this->n_patches() == 0) {
    return;
  }

  swap_alt();
  resize(this->n_prts);

  int n_blocks = (this->n_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  if (n_blocks > 32768)
    n_blocks = 32768;
  dim3 dimGrid(n_blocks);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_reorder_and_offsets<BS><<<dimGrid, dimBlock>>>(
    *this, this->n_prts, d_idx.data().get(), d_id.data().get(),
    d_off.data().get(), this->n_blocks);
  cuda_sync_if_enabled();

  need_reorder = false;
}

// ----------------------------------------------------------------------
// k_reorder

template <typename BS>
__global__ static void k_reorder(DMparticlesCuda<BS> dmprts, int n_prts,
                                 const uint* d_ids)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i < n_prts) {
    dmprts.storage.store(dmprts.alt_storage[d_ids[i]],
                         i); // storage[i] = alt_storage[j]
  }
}

// ----------------------------------------------------------------------
// reorder

template <typename BS>
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

template <typename BS>
void cuda_mparticles<BS>::reorder(const psc::device_vector<uint>& d_id)
{
  if (this->n_prts == 0) {
    return;
  }

  swap_alt();
  resize(this->n_prts);

  dim3 dimGrid((this->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);

  k_reorder<BS>
    <<<dimGrid, THREADS_PER_BLOCK>>>(*this, this->n_prts, d_id.data().get());
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// setup_internals

template <typename BS>
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
// size

template <typename BS>
uint cuda_mparticles<BS>::size()
{
  return this->n_prts;
}

// ----------------------------------------------------------------------
// inject_initial
//
// adds particles initially, ie., into an empty cmprts
// does not complete setting correct internal state
// (setup_internal() needs to be called next)

template <typename BS>
void cuda_mparticles<BS>::inject_initial(
  const std::vector<Particle>& buf, const std::vector<uint>& n_prts_by_patch)
{
  thrust::host_vector<uint> h_off(this->by_block_.d_off);

  assert(this->storage.size() == 0);
  assert(this->n_prts == 0);

  uint buf_n = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    assert(h_off[p * this->n_blocks_per_patch] == 0);
    assert(h_off[(p + 1) * this->n_blocks_per_patch] == 0);
    buf_n += n_prts_by_patch[p];
  }

  resize(buf_n);

  HMparticlesCudaStorage h_storage{buf_n};

  auto it = buf.begin();
  uint off = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    auto n_prts = n_prts_by_patch[p];
    h_off[p * this->n_blocks_per_patch] = off;
    h_off[(p + 1) * this->n_blocks_per_patch] = off + n_prts;

    for (int n = 0; n < n_prts; n++) {
      auto prt = *it++;
      this->checkInPatchMod(prt.x);
      h_storage.store(prt, off + n);
    }

    off += n_prts;
  }
  this->n_prts = off;

  this->storage = h_storage;
  thrust::copy(h_off.begin(), h_off.end(), this->by_block_.d_off.begin());
}

// ----------------------------------------------------------------------
// inject

template <typename BS>
void cuda_mparticles<BS>::inject(const std::vector<Particle>& buf,
                                 const std::vector<uint>& buf_n_by_patch)
{
  if (this->n_prts == 0) {
    // if there are no particles yet, we basically just initialize from the
    // buffer
    inject_initial(buf, buf_n_by_patch);
    setup_internals();
    return;
  }

  using Double3 = Vec3<double>;

  uint buf_n = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    buf_n += buf_n_by_patch[p];
    //    printf("p %d buf_n_by_patch %d\n", p, buf_n_by_patch[p]);
  }
  //  printf("buf_n %d\n", buf_n);

  HMparticlesCudaStorage h_storage(buf_n);
  thrust::host_vector<uint> h_bidx(buf_n);
  // thrust::host_vector<uint> h_id(buf_n);

  uint off = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    for (int n = 0; n < buf_n_by_patch[p]; n++) {
      auto prt = buf[off + n];
      h_storage.store(prt, off + n);
      auto bidx = this->blockIndex(prt, p);
      assert(bidx >= 0 && bidx < this->n_blocks);
      h_bidx[off + n] = bidx;
      ;
      // h_id[off + n] = this->n_prts + off + n;
    }
    off += buf_n_by_patch[p];
  }
  assert(off == buf_n);

  if (need_reorder) {
    reorder();
  }

  // assert(check_in_patch_unordered_slow());

  this->by_block_.find_indices_ids(*this);
  // assert(check_bidx_id_unordered_slow());

  resize(this->n_prts + buf_n);

  thrust::copy(h_storage.begin(), h_storage.end(),
               this->storage.begin() + this->n_prts);
  thrust::copy(h_bidx.begin(), h_bidx.end(),
               this->by_block_.d_idx.begin() + this->n_prts);
  // thrust::copy(h_id.begin(), h_id.end(), d_id + n_prts);
  // FIXME, looks like ids up until n_prts have already been set above
  thrust::sequence(this->by_block_.d_id.data(),
                   this->by_block_.d_id.data() + this->n_prts + buf_n);

  // assert(check_ordered());

  this->n_prts += buf_n;

  this->by_block_.stable_sort();

  this->by_block_.reorder_and_offsets(*this);

  // assert(check_ordered());
}

// ----------------------------------------------------------------------
// get_offsets

template <typename BS>
std::vector<uint> cuda_mparticles<BS>::get_offsets() const
{
  thrust::host_vector<uint> h_off(this->by_block_.d_off);
  std::vector<uint> off(this->n_patches() + 1);
  for (int p = 0; p <= this->n_patches(); p++) {
    off[p] = h_off[p * this->n_blocks_per_patch];
  }
  return off;
}

inline void copy(const MparticlesCudaStorage& from, HMparticlesCudaStorage& to)
{
  assert(from.size() == to.size());
  thrust::copy(from.xi4.begin(), from.xi4.end(), to.xi4.begin());
  thrust::copy(from.pxi4.begin(), from.pxi4.end(), to.pxi4.begin());
}

// ----------------------------------------------------------------------
// get_particles

template <typename BS>
std::vector<typename cuda_mparticles<BS>::Particle>
cuda_mparticles<BS>::get_particles()
{
  reorder(); // FIXME? by means of this, this function disturbs the state...

  assert(this->storage.size() == this->n_prts);
  HMparticlesCudaStorage h_storage{this->n_prts};
  copy(this->storage, h_storage);

  std::vector<Particle> prts;
  prts.reserve(this->n_prts);

  for (int n = 0; n < this->n_prts; n++) {
    prts.emplace_back(h_storage[n]);
  }

  return prts;
}

#include "cuda_mparticles_gold.cxx"
#include "cuda_mparticles_checks.cxx"

template struct cuda_mparticles<BS144>;
template struct cuda_mparticles<BS444>;
