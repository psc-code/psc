
#include "cuda_mparticles.hxx"
#include "cuda_bits.h"

// ----------------------------------------------------------------------
// check_in_patch_unordered_slow

template <typename BS>
bool cuda_mparticles<BS>::check_in_patch_unordered_slow()
{
  auto n_prts_by_patch = this->sizeByPatch();

  uint off = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = this->blockIndex(this->storage[off + n], p);
      if (!(bidx >= 0 && bidx <= this->n_blocks))
        return false;
    }
    off += n_prts_by_patch[p];
  }

  if (!(off == this->n_prts))
    return false;
  // printf("PASS: cuda_mparticles_check_in_patch_unordered_slow()\n");
  return true;
}

// ----------------------------------------------------------------------
// check_bix_id_unordered_slow
//
// checks that block indices are correct,
// id is just enumerating particles

template <typename BS>
bool cuda_mparticles<BS>::check_bidx_id_unordered_slow()
{
  auto n_prts_by_patch = this->sizeByPatch();

  uint off = 0;
  for (int p = 0; p < this->n_patches(); p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      int bidx = this->blockIndex(this->storage[off + n], p);
      if (!(bidx == this->by_block_.d_idx[off + n]))
        return false;
      if (!(off + n == this->by_block_.d_id[off + n]))
        return false;
    }
    off += n_prts_by_patch[p];
  }

  if (!(off == this->n_prts))
    return false;
  // printf("PASS: cuda_mparticles_check_bidx_id_unordered_slow()\n");
  return true;
}

// ----------------------------------------------------------------------
// check_ordered

template <typename BS>
bool cuda_mparticles<BS>::check_ordered()
{
  HMparticlesCudaStorage h_storage(this->storage);
  thrust::host_vector<uint> h_off(this->by_block_.d_off);
  thrust::host_vector<uint> h_id(this->by_block_.d_id);

  // printf("check_ordered: need_reorder %s\n", need_reorder ? "true" :
  // "false");

  uint off = 0;
  for (int b = 0; b < this->n_blocks; b++) {
    int p = b / this->n_blocks_per_patch;
    uint off_b = h_off[b], off_e = h_off[b + 1];
    if (!(off_e >= off_b))
      return false;
    // printf("check_ordered: block %d: %d -> %d (patch %d)\n", b, off_b, off_e,
    // p);
    if (!(off_b == off))
      return false;
    for (int n = h_off[b]; n < h_off[b + 1]; n++) {
      DParticleCuda prt;
      if (need_reorder) {
        prt = h_storage[h_id[n]];
      } else {
        prt = h_storage[n];
      }
      uint bidx = this->blockIndex(prt, p);
      // printf("check_ordered: bidx %d\n", bidx);
      if (b != bidx) {
        printf("check_ordered: b %d bidx %d n %d p %d x %g %g %g\n", b, bidx, n,
               p, prt.x[0], prt.x[1], prt.x[2]);
        Int3 bpos = this->blockPosition(prt.x);
        printf("block_pos %d %d\n", bpos[1], bpos[2]);
      }
      if (!(b == bidx))
        return false;
    }
    off += off_e - off_b;
  }
  if (!(off == this->n_prts))
    return false;
  // printf("PASS: cuda_mparticles_check_ordered:\n");
  return true;
}

// ----------------------------------------------------------------------
// check_bidx_after_push

template <typename BS>
bool cuda_mparticles<BS>::check_bidx_after_push()
{
  bool ok = true;

  thrust::host_vector<uint> h_off(this->by_block_.d_off);
  thrust::host_vector<uint> h_bidx(this->by_block_.d_idx);
  HMparticlesCudaStorage h_storage(this->storage);

  for (int p = 0; p < this->n_patches(); p++) {
    int begin = h_off[p * this->n_blocks_per_patch];
    int end = h_off[(p + 1) * this->n_blocks_per_patch];
    for (int n = begin; n < end; n++) {
      DParticleCuda prt = h_storage[n];
      int bidx = h_bidx[n];
      int bidx2 = this->blockIndex(prt, p);
      if (bidx2 < 0)
        bidx2 = this->n_blocks;
      if (bidx != bidx2) {
        mprintf("check_bidx: n %d: x %g %g %g bidx %d/%d\n", n, prt.x[0],
                prt.x[1], prt.x[2], bidx, bidx2);
        ok = false;
      }
    }
  }
  return ok;
}
