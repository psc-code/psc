
#include "cuda_mparticles.h"

// ======================================================================
// cuda_mparticles_base

// ----------------------------------------------------------------------
// ctor

cuda_mparticles_base::cuda_mparticles_base(const Grid_t& grid, const Int3& bs)
  : grid_(grid),
  n_patches(grid.patches.size())
{
  this->bs = bs;
  
  for (int d = 0; d < 3; d++) {
    assert(grid.ldims[d] % bs[d] == 0);
    b_mx[d] = grid.ldims[d] / bs[d];
  }
  
  n_blocks_per_patch = b_mx[0] * b_mx[1] * b_mx[2];
  n_blocks = n_patches * n_blocks_per_patch;
  
  d_off.resize(n_blocks + 1);
}

// ----------------------------------------------------------------------
// reserve_all
//
// FIXME, eventually should do its own thing, but for now we better keep the size
// the same as for the rest of the arrays
  
void cuda_mparticles_base::reserve_all()
{
  d_xi4.resize(n_alloced);
  d_pxi4.resize(n_alloced);
}

// ----------------------------------------------------------------------
// resize_all
//
// FIXME, this function currently is used in two contexts:
// - to implement mprts::resize_all(), but in this case we
//   need to be careful. It's destructive, which is unexpected.
//   we might want to only support (and check for) the case of
//   resizing from 0 size.
//   in this case, we also should check that things fit into what's
//   alloced (also: a very similar issues is cuda_mparticles_reserve_all()
//   which doesn't realloc but destroy, again that's unexpected behavior
// - to reset the internal n_prts_by_patch as part of sorting etc.
//   in that case, we supposedly know what we're doing, so we at most need
//   to check that we aren't beyond our allocated space
  
void cuda_mparticles_base::resize_all(const uint *n_prts_by_patch)
{
  thrust::host_vector<uint> h_off(this->n_blocks + 1);
    
  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    h_off[p * n_blocks_per_patch] = off;
    off += n_prts_by_patch[p];
    // printf("set_n_prts p%d: %d\n", p, n_prts_by_patch[p]);
  }
  h_off[n_blocks] = off;
  n_prts = off;
    
  thrust::copy(h_off.begin(), h_off.end(), d_off.begin());
}

