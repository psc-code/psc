
#pragma once

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
  
  DParticleCuda load(int n) const
  {
    float4 _xi4 = xi4[n];
    float4 _pxi4 = pxi4[n];
    return {{_xi4.x, _xi4.y, _xi4.z}, {_pxi4.x, _pxi4.y, _pxi4.z}, _pxi4.w, cuda_float_as_int(_xi4.w)};
  }

  void store(const DParticleCuda& prt, int n)
  {
    float4 _xi4 = { prt.x()[0], prt.x()[1], prt.x()[2], cuda_int_as_float(prt.kind()) };
    float4 _pxi4 = { prt.u()[0], prt.u()[1], prt.u()[2], prt.qni_wni() };
    xi4[n] = _xi4;
    pxi4[n] = _pxi4;
  }
  
  thrust::device_vector<float4> xi4;
  thrust::device_vector<float4> pxi4;
};

// ======================================================================
// HMparticlesCudaStorage

struct HMparticlesCudaStorage
{
  HMparticlesCudaStorage(size_t n)
    : xi4{n}, pxi4{n}
  {}

  HMparticlesCudaStorage(const MparticlesCudaStorage& storage)
    : xi4{storage.xi4}, pxi4{storage.pxi4}
  {}

  // FIXME, why so many warnings?
  void resize(size_t n)
  {
    xi4.resize(n);
    pxi4.resize(n);
  }
  
  DParticleCuda load(int n) const
  {
    float4 _xi4 = xi4[n];
    float4 _pxi4 = pxi4[n];
    return {{_xi4.x, _xi4.y, _xi4.z}, {_pxi4.x, _pxi4.y, _pxi4.z}, _pxi4.w, cuda_float_as_int(_xi4.w)};
  }

  void store(const DParticleCuda& prt, int n)
  {
    float4 _xi4 = { prt.x()[0], prt.x()[1], prt.x()[2], cuda_int_as_float(prt.kind()) };
    float4 _pxi4 = { prt.u()[0], prt.u()[1], prt.u()[2], prt.qni_wni() };
    xi4[n] = _xi4;
    pxi4[n] = _pxi4;
  }
  
  thrust::host_vector<float4> xi4;
  thrust::host_vector<float4> pxi4;
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
  void store(const DParticleCuda& prt, int n)
  {
    store_position(prt, n);
    store_momentum(prt, n);
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

protected:
  void resize(uint size);

public:
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
  using particle_t = DParticleCuda;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;
  using DMparticles = DMparticlesCuda<BS>;

  struct Patch
  {
    struct injector
    {
      injector(const Patch& patch)
	: patch_{patch},
	  n_prts_{0}
      {
	assert(patch_.p_ == patch_.cmprts_.injector_n_prts_by_patch_.size());
      }

      ~injector()
      {
	auto& cmprts = patch_.cmprts_;
	cmprts.injector_n_prts_by_patch_.push_back(n_prts_);
	if (patch_.p_ == cmprts.n_patches - 1) {
	  cmprts.inject(cmprts.injector_buf_, cmprts.injector_n_prts_by_patch_);
	  cmprts.injector_n_prts_by_patch_.clear();
	  cmprts.injector_buf_.clear();
	}
      }
      
      void raw(const particle_t& prt)
      {
	patch_.cmprts_.injector_buf_.push_back(prt);
	n_prts_++;
      }

      void operator()(const particle_inject& new_prt)
      {
	using Double3 = Vec3<double>;
	
	auto& cmprts = patch_.cmprts_;
	auto& patch = cmprts.grid_.patches[patch_.p_];
	auto x = Double3::fromPointer(new_prt.x) - patch.xb;
	auto prt = particle_t{Real3(x), Real3(Double3::fromPointer(new_prt.u)),
			      real_t(new_prt.w), new_prt.kind};
	patch_.cmprts_.injector_buf_.push_back(prt);
	n_prts_++;
      }

      void operator()(const std::vector<particle_t>& buf)
      {
	auto& injector_buf = patch_.cmprts_.injector_buf_;
	injector_buf.insert(injector_buf.end(), buf.begin(), buf.end());
	n_prts_ += buf.size();
      }
      
    private:
      const Patch patch_;
      uint n_prts_;
    };

    Patch(cuda_mparticles& cmprts, int p)
      : cmprts_(cmprts), p_(p)
    {}

    injector injector() { return {*this}; }

  private:
    cuda_mparticles& cmprts_;
    int p_;
  };

  cuda_mparticles(const Grid_t& grid);

  Patch operator[](int p) { return Patch{*this, p}; }
  

  uint get_n_prts();
  void inject(const std::vector<particle_t>& buf, const std::vector<uint>& buf_n_by_patch);

  std::vector<particle_t> get_particles(int beg, int end);
  std::vector<particle_t> get_particles(int p);

  uint start(int p);

  void dump(const std::string& filename) const;
  void dump_by_patch(uint *n_prts_by_patch);

  // internal / testing use
  void inject_initial(const std::vector<particle_t>& buf,
		      const std::vector<uint>& n_prts_by_patch);
  void setup_internals();

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

private:
  std::vector<particle_t> injector_buf_;
  std::vector<uint> injector_n_prts_by_patch_;
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

