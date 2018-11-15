
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"
#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles_sort.cuh"

#include "particles.hxx"
#include "psc_bits.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// ======================================================================
// MparticlesCudaStorage

struct MparticlesCudaStorage
{
  void resize(size_t n)
  {
    xi4.resize(n);
    pxi4.resize(n);
  }
  
  thrust::device_vector<float4> xi4;
  thrust::device_vector<float4> pxi4;
};

// ======================================================================
// DMparticlesCudaStorage

struct DMparticlesCudaStorage
{
  // FIXME, could be operator[]
  __host__ __device__
  DParticleCuda load(int n)
  {
    float4 _xi4 = xi4[n];
    float4 _pxi4 = pxi4[n];
    return {{_xi4.x, _xi4.y, _xi4.z}, {_pxi4.x, _pxi4.y, _pxi4.z}, _pxi4.w, cuda_float_as_int(_xi4.w)};
  }

  __host__ __device__
  void store_position(const DParticleCuda& prt, int n)
  {
    float4 _xi4 = { prt.x()[0], prt.x()[1], prt.x()[2], cuda_int_as_float(prt.kind()) };
    xi4[n] = _xi4;
  }

  __host__ __device__
  void store_momentum(const DParticleCuda& prt, int n)
  {
    float4 _pxi4 = { prt.u()[0], prt.u()[1], prt.u()[2], prt.qni_wni() };
    pxi4[n] = _pxi4;
  }
  
  float4 *xi4;
  float4 *pxi4;
};

// ======================================================================
// cuda_mparticles_base

template<typename BS>
struct cuda_mparticles_base : cuda_mparticles_indexer<BS>
{
  cuda_mparticles_base(const Grid_t& grid);
  // copy constructor would work fine, but we don't want to copy everything
  // by accident
  cuda_mparticles_base(const cuda_mparticles<BS>&) = delete;

  void reserve_all(uint sinze);
  void resize_all(const uint *n_prts_by_patch);
  void get_size_all(uint *n_prts_by_patch);

  // per particle
  MparticlesCudaStorage storage;

  // per block
  cuda_mparticles_sort2 by_block_;

  uint n_prts = 0;                       // total # of particles across all patches
  const Grid_t& grid_;
};

// ----------------------------------------------------------------------
// cuda_mparticles

template<typename BS>
struct DMparticlesCuda;

template<typename _BS>
struct cuda_mparticles : cuda_mparticles_base<_BS>
{
  using BS = _BS;
  using particle_t = BndpParticleCuda;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;
  using DMparticles = DMparticlesCuda<BS>;

  cuda_mparticles(const Grid_t& grid);

  void reserve_all(const uint *n_prts_by_patch);
  uint get_n_prts();
  void setup_internals();
  void inject_buf(const cuda_mparticles_prt *buf, const uint *buf_n_by_patch);
  void inject_buf(const particle_inject *buf, const uint *buf_n_by_patch);

  std::vector<cuda_mparticles_prt> get_particles(int beg, int end);
  std::vector<cuda_mparticles_prt> get_particles(int p);

  uint start(int p);
  
  template<typename F>
  void set_particles(uint p, F getter);

  void dump(const std::string& filename) const;
  void dump_by_patch(uint *n_prts_by_patch);

public:
  void find_block_indices_ids(thrust::device_vector<uint>& d_idx, thrust::device_vector<uint>& d_id);
  void find_cell_indices_ids(thrust::device_vector<uint>& d_idx, thrust::device_vector<uint>& d_id);
  void reorder();
  void reorder(const thrust::device_vector<uint>& d_id);
  void reorder_and_offsets(const thrust::device_vector<uint>& d_idx, const thrust::device_vector<uint>& d_id,
			   thrust::device_vector<uint>& d_off);
  void reorder_and_offsets_slow();
  void swap_alt();

  bool check_in_patch_unordered_slow();
  bool check_bidx_id_unordered_slow();
  bool check_ordered();
  bool check_bidx_after_push();

  void resize(uint n_prts);

public:
  MparticlesCudaStorage alt_storage; // storage for out-of-place reordering of particle data

  std::vector<Real3> xb_by_patch; // lower left corner for each patch

  bool need_reorder = { false };            // particles haven't yet been put into their sorted order
};

template<typename BS_>
struct DMparticlesCuda : DParticleIndexer<BS_>
{
  using BS = BS_;
  using typename DParticleIndexer<BS>::real_t;
  
  static const int MAX_N_KINDS = 4;

  DMparticlesCuda(cuda_mparticles<BS>& cmprts)
    : DParticleIndexer<BS>{cmprts},
      dt_(cmprts.grid_.dt),
      fnqs_(cmprts.grid_.norm.fnqs),
      fnqxs_(cmprts.grid_.domain.dx[0] * fnqs_ / dt_),
      fnqys_(cmprts.grid_.domain.dx[1] * fnqs_ / dt_),
      fnqzs_(cmprts.grid_.domain.dx[2] * fnqs_ / dt_),
      dqs_(.5f * cmprts.grid_.norm.eta * dt_),
      storage{cmprts.storage.xi4.data().get(), cmprts.storage.pxi4.data().get()},
      alt_storage{cmprts.alt_storage.xi4.data().get(), cmprts.alt_storage.pxi4.data().get()},
      off_(cmprts.by_block_.d_off.data().get()),
      bidx_(cmprts.by_block_.d_idx.data().get()),
      id_(cmprts.by_block_.d_id.data().get()),
      n_blocks_(cmprts.n_blocks)
  {
    auto& grid = cmprts.grid_;
    
    int n_kinds = grid.kinds.size();
    assert(n_kinds <= MAX_N_KINDS);
    for (int k = 0; k < n_kinds; k++) {
      dq_[k] = dqs_ * grid.kinds[k].q / grid.kinds[k].m;
      q_inv_[k] = 1.f / grid.kinds[k].q;
      q_[k] = grid.kinds[k].q;
      m_[k] = grid.kinds[k].m;
    }
  }

  __device__ real_t dt() const { return dt_; }
  __device__ real_t fnqs() const { return fnqs_; }
  __device__ real_t fnqxs() const { return fnqxs_; }
  __device__ real_t fnqys() const { return fnqys_; }
  __device__ real_t fnqzs() const { return fnqzs_; }
  __device__ real_t q_inv(int k) const { return q_inv_[k]; }
  __device__ real_t dq(int k) const { return dq_[k]; }
  __device__ real_t q(int k) const { return q_[k]; }
  __device__ real_t m(int k) const { return m_[k]; }

private:
  real_t dt_;
  real_t fnqs_;
  real_t fnqxs_, fnqys_, fnqzs_;
  real_t dqs_;
  real_t dq_[MAX_N_KINDS];
  real_t q_inv_[MAX_N_KINDS];
  real_t q_[MAX_N_KINDS];
  real_t m_[MAX_N_KINDS];
public:
  DMparticlesCudaStorage storage;
  DMparticlesCudaStorage alt_storage;
  uint *off_;
  uint *bidx_;
  uint *id_;
  uint n_blocks_;
};

// ======================================================================
// cuda_mparticles implementation

// ----------------------------------------------------------------------
// set_particles

template<typename BS>
template<typename F>
void cuda_mparticles<BS>::set_particles(uint p, F getter)
{
  // FIXME, doing the copy here all the time would be nice to avoid
  // making sue we actually have a valid d_off would't hurt, either
  thrust::host_vector<uint> h_off(this->by_block_.d_off);

  uint off = h_off[p * this->n_blocks_per_patch];
  uint n_prts = h_off[(p+1) * this->n_blocks_per_patch] - off;
  
  thrust::host_vector<float4> xi4(n_prts);
  thrust::host_vector<float4> pxi4(n_prts);

  for (int n = 0; n < n_prts; n++) {
    cuda_mparticles_prt prt = getter(n);
    this->checkInPatchMod(prt.x());

    xi4[n].x  = prt.x()[0];
    xi4[n].y  = prt.x()[1];
    xi4[n].z  = prt.x()[2];
    xi4[n].w  = cuda_int_as_float(prt.kind());
    pxi4[n].x = prt.u()[0];
    pxi4[n].y = prt.u()[1];
    pxi4[n].z = prt.u()[2];
    pxi4[n].w = prt.qni_wni();
  }

  thrust::copy(xi4.begin(), xi4.end(), &this->storage.xi4[off]);
  thrust::copy(pxi4.begin(), pxi4.end(), &this->storage.pxi4[off]);
}



#endif
