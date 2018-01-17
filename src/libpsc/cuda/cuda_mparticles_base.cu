
#include "cuda_mparticles.h"

#include "psc_bits.h"
#include "cuda_bits.h"

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
    b_dxi[d] = 1.f / (bs[d] * grid.dx[d]);
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

// ----------------------------------------------------------------------
// get_size_all

void cuda_mparticles_base::get_size_all(uint *n_prts_by_patch)
{
  thrust::host_vector<uint> h_off(d_off);

  for (int p = 0; p < n_patches; p++) {
    n_prts_by_patch[p] = h_off[(p+1) * n_blocks_per_patch] - h_off[p * n_blocks_per_patch];
    //printf("p %d n_prts_by_patch %d\n", p, n_prts_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// to_device

void cuda_mparticles_base::to_device(float4 *xi4, float4 *pxi4,
				     uint n_prts, uint off)
{
  assert(off + n_prts <= n_alloced);
  thrust::copy(xi4, xi4 + n_prts, d_xi4.begin() + off);
  thrust::copy(pxi4, pxi4 + n_prts, d_pxi4.begin() + off);
}

// ----------------------------------------------------------------------
// from_device

void cuda_mparticles_base::from_device(float4 *xi4, float4 *pxi4,
				       uint n_prts, uint off)
{
  assert(off + n_prts <= n_alloced);
  thrust::copy(d_xi4.begin() + off, d_xi4.begin() + off + n_prts, xi4);
  thrust::copy(d_pxi4.begin() + off, d_pxi4.begin() + off + n_prts, pxi4);
}

// ----------------------------------------------------------------------
// set_particles

void cuda_mparticles_base::set_particles(uint n_prts, uint off,
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

