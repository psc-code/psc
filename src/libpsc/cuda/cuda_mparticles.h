
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"

#include "grid.hxx"
#include "particles.hxx"
#include "psc_bits.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>

// ======================================================================
// cuda_mparticles_indexer

struct cuda_mparticles_indexer
{
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_indexer() = default; // FIXME, should go away
  
  cuda_mparticles_indexer(const Int3& b_mx, const Real3& b_dxi)
    : b_mx_(b_mx), b_dxi_(b_dxi)
  {}
  
  Int3 b_mx_;          // number of blocks per direction in each patch
  Real3 b_dxi_;        // inverse of block size (in actual length units)
};

// ======================================================================
// cuda_mparticles_base

struct cuda_mparticles_base
{
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_base(const Grid_t& grid, const Int3& bs);
  // copy constructor would work fine, be don't want to copy everything
  // by accident
  cuda_mparticles_base(const cuda_mparticles&) = delete;

  void reserve_all();
  void resize_all(const uint *n_prts_by_patch);
  void get_size_all(uint *n_prts_by_patch);

  // protected:
  void to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  
  // per particle
  thrust::device_vector<float4> d_xi4;
  thrust::device_vector<float4> d_pxi4;

  // per block
  thrust::device_vector<uint> d_off;     // particles per block
                                         // are at indices [offsets[block] .. offsets[block+1]-1[

  uint n_prts = {};                      // total # of particles across all patches
  uint n_alloced = {};                   // size of particle-related arrays as allocated
  Int3 bs;                               // block size
  uint n_blocks_per_patch;               // number of blocks per patch
  uint n_blocks;                         // number of blocks in all patches in mprts
  cuda_mparticles_indexer indexer;

  const uint n_patches;
  const Grid_t& grid_;
};

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles : cuda_mparticles_base
{
  cuda_mparticles(const Grid_t& grid, const Int3& bs);

  void reserve_all(const uint *n_prts_by_patch);
  uint get_n_prts();
  void setup_internals();
  void inject(const cuda_mparticles_prt *buf, uint *buf_n_by_patch);
  const particle_cuda_real_t *patch_get_b_dxi(int p);
  const int *patch_get_b_mx(int p);

  template<typename F>
  void set_particles(uint p, F getter);

  template<typename F>
  void get_particles(uint p, F setter);
  
  void dump();
  void dump_by_patch(uint *n_prts_by_patch);

public:
  void find_block_indices_ids();
  void stable_sort_by_key();
  void reorder();
  void reorder_and_offsets();
  void reorder_and_offsets_slow();
  int get_block_idx(float4 xi4, int p);
  void swap_alt();

  bool check_in_patch_unordered_slow();
  bool check_bidx_id_unordered_slow();
  bool check_ordered();
  
public:
  // per particle
  thrust::device_vector<float4> d_alt_xi4;  // storage for out-of-place reordering of particle data
  thrust::device_vector<float4> d_alt_pxi4;
  thrust::device_vector<uint> d_bidx;       // block index (incl patch) per particle
  thrust::device_vector<uint> d_id;         // particle id for sorting

  std::vector<Real3> xb_by_patch; // lower left corner for each patch

  bool need_reorder = { false };            // particles haven't yet been put into their sorted order
};

// ======================================================================
// cuda_mparticles implementation

// ----------------------------------------------------------------------
// set_particles

template<typename F>
void cuda_mparticles::set_particles(uint p, F getter)
{
  // FIXME, doing the copy here all the time would be nice to avoid
  // making sue we actually have a valid d_off would't hurt, either
  thrust::host_vector<uint> h_off(d_off);

  uint off = h_off[p * n_blocks_per_patch];
  uint n_prts = h_off[(p+1) * n_blocks_per_patch] - off;
  
  thrust::host_vector<float4> xi4(n_prts);
  thrust::host_vector<float4> pxi4(n_prts);

  for (int n = 0; n < n_prts; n++) {
    struct cuda_mparticles_prt prt = getter(n);

    for (int d = 0; d < 3; d++) {
      int bi = fint(prt.xi[d] * indexer.b_dxi_[d]);
      if (bi < 0 || bi >= indexer.b_mx_[d]) {
	printf("XXX xi %g %g %g\n", prt.xi[0], prt.xi[1], prt.xi[2]);
	printf("XXX n %d d %d xi4[n] %g biy %d // %d\n",
	       n, d, prt.xi[d], bi, indexer.b_mx_[d]);
	if (bi < 0) {
	  prt.xi[d] = 0.f;
	} else {
	  prt.xi[d] *= (1. - 1e-6);
	}
      }
      bi = floorf(prt.xi[d] * indexer.b_dxi_[d]);
      assert(bi >= 0 && bi < indexer.b_mx_[d]);
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

  thrust::copy(xi4.begin(), xi4.end(), &d_xi4[off]);
  thrust::copy(pxi4.begin(), pxi4.end(), &d_pxi4[off]);
}

// ----------------------------------------------------------------------
// get_particles

template<typename F>
void cuda_mparticles::get_particles(uint p, F setter)
{
  // FIXME, doing the copy here all the time would be nice to avoid
  // making sue we actually have a valid d_off would't hurt, either
  thrust::host_vector<uint> h_off(d_off);

  uint off = h_off[p * n_blocks_per_patch];
  uint n_prts = h_off[(p+1) * n_blocks_per_patch] - off;

  reorder();

  thrust::host_vector<float4> xi4(&d_xi4[off], &d_xi4[off + n_prts]);
  thrust::host_vector<float4> pxi4(&d_pxi4[off], &d_pxi4[off + n_prts]);

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

    setter(n, prt);

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

}



#endif
