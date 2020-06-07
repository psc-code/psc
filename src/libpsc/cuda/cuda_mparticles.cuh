
#pragma once

#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles_sort.cuh"
#include "injector_buffered.hxx"
#include "mparticles_patch_cuda.hxx"

#include "particles.hxx"
#include "psc_bits.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// ======================================================================
// DParticleProxy

template<typename DMparticlesCuda>
struct DParticleProxy
{
  using real_t = DParticleCuda::real_t;
  using Real3 = DParticleCuda::Real3;
  
  __device__
  DParticleProxy(const DParticleCuda& prt, const DMparticlesCuda& dmprts)
    : prt_{prt}, dmprts_{dmprts}
  {}

  __device__ Real3  x() const { return prt_.x; }
  __device__ Real3& x()       { return prt_.x; }
  __device__ Real3  u() const { return prt_.u; }
  __device__ Real3& u()       { return prt_.u; }
  __device__ int kind() const { return prt_.kind; }
  __device__ real_t qni_wni() const { return prt_.qni_wni; }

  __device__ real_t q() const { return dmprts_.q(prt_.kind); }
  __device__ real_t m() const { return dmprts_.m(prt_.kind); }
  
private:
  DParticleCuda prt_;
  const DMparticlesCuda& dmprts_;
};

// ======================================================================
// ParticleCudaStorage

struct ParticleCudaStorage
{
  __host__ __device__
  ParticleCudaStorage(float4 xi4, float4 pxi4)
    : xi4{xi4}, pxi4{pxi4}
  {}

  __host__ __device__
  ParticleCudaStorage(const DParticleCuda& prt)
  : xi4{prt.x[0], prt.x[1], prt.x[2], cuda_int_as_float(prt.kind)},
    pxi4{prt.u[0], prt.u[1], prt.u[2], prt.qni_wni}
  {}

  template<typename DParticleProxy>
  __device__
  ParticleCudaStorage(const DParticleProxy& prt)
  : xi4{prt.x()[0], prt.x()[1], prt.x()[2], cuda_int_as_float(prt.kind())},
    pxi4{prt.u()[0], prt.u()[1], prt.u()[2], prt.qni_wni()}
  {}

  __host__ __device__
  operator DParticleCuda()
  {
    return {{xi4.x, xi4.y, xi4.z}, {pxi4.x, pxi4.y, pxi4.z}, pxi4.w, cuda_float_as_int(xi4.w), psc::particle::Id(), psc::particle::Tag()};
  }
  
  float4 xi4;
  float4 pxi4;
};

// ======================================================================
// MparticlesCudaStorage_

template<typename T>
struct MparticlesCudaStorage_
{
  MparticlesCudaStorage_() = default;

  __host__
  MparticlesCudaStorage_(const T& xi4, const T& pxi4)
    : xi4{xi4}, pxi4{pxi4}
  {}
  
  template<typename OtherStorage>
  __host__
  MparticlesCudaStorage_(OtherStorage& other)
  : xi4{other.xi4}, pxi4{other.pxi4}
  {}

  __host__
  MparticlesCudaStorage_(uint n)
  : xi4(n), pxi4(n)
  {}

  __host__
  void resize(size_t n)
  {
    xi4.resize(n);
    pxi4.resize(n);
  }
  
  // FIXME, could be operator[]

  __device__
  DParticleCuda load_device(int n) const { return ParticleCudaStorage{xi4[n], pxi4[n]}; }

  __host__
  DParticleCuda load(int n) const { return ParticleCudaStorage{xi4[n], pxi4[n]}; }

  __host__ __device__
  void store(const DParticleCuda& prt, int n)
  {
    auto st = ParticleCudaStorage{prt};
    xi4[n] = st.xi4;
    pxi4[n] = st.pxi4;
  }

  T xi4;
  T pxi4;
};

// ======================================================================
// MparticlesCudaStorage

using MparticlesCudaStorage = MparticlesCudaStorage_<thrust::device_vector<float4>>;

// ======================================================================
// HMparticlesCudaStorage

using HMparticlesCudaStorage = MparticlesCudaStorage_<thrust::host_vector<float4>>;

// ======================================================================
// DMparticlesCudaStorage

struct DMparticlesCudaStorage : MparticlesCudaStorage_<float4*>
{
  using Base = MparticlesCudaStorage_<float4*>;
  using Base::Base;

  template<typename DMparticlesCuda>
  __device__
  void store_position(const DParticleProxy<DMparticlesCuda>& prt, int n)
  {
    auto st = ParticleCudaStorage{prt};
    xi4[n] = st.xi4;
  }

  template<typename DMparticlesCuda>
  __device__
  void store_momentum(const DParticleProxy<DMparticlesCuda>& prt, int n)
  {
    auto st = ParticleCudaStorage{prt};
    pxi4[n] = st.pxi4;
  }

  template<typename DMparticlesCuda>
  __device__
  DParticleProxy<DMparticlesCuda> load_proxy(const DMparticlesCuda& dmprts, int n)
  {
    return {load_device(n), dmprts};
  }
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
  std::vector<uint> sizeByPatch() const;
  void clear() { storage.resize(0); by_block_.clear(); n_prts = 0; }

  // per particle
  MparticlesCudaStorage storage;

  // per block
  cuda_mparticles_sort_by_block by_block_;

  uint n_prts = 0;                       // total # of particles across all patches
  const Grid_t& grid_;
};

// ======================================================================
// cuda_mparticles

template<typename BS>
struct DMparticlesCuda;

template<typename _BS>
struct cuda_mparticles : cuda_mparticles_base<_BS>
{
  using BS = _BS;
  using Particle = DParticleCuda;
  using real_t = Particle::real_t;
  using Real3 = Vec3<real_t>;
  using DMparticles = DMparticlesCuda<BS>;
  using Patch = ConstPatchCuda<cuda_mparticles>;
  using BndpParticle = DParticleCuda;
  using BndBuffer = std::vector<BndpParticle>;
  using BndBuffers = std::vector<BndBuffer>;

  cuda_mparticles(const Grid_t& grid);

  InjectorBuffered<cuda_mparticles> injector() { return {*this}; }
  ConstAccessorCuda_<cuda_mparticles> accessor() { return {*this}; }
  
  uint size();
  void inject(const std::vector<Particle>& buf, const std::vector<uint>& buf_n_by_patch);

  std::vector<uint> get_offsets() const;
  std::vector<Particle> get_particles();
  std::vector<Particle> get_particles(int p);
  Particle get_particle(int p, int n);

  void dump(const std::string& filename) const;
  void dump_by_patch(uint *n_prts_by_patch);

  // internal / testing use
  void inject_initial(const std::vector<Particle>& buf,
		      const std::vector<uint>& n_prts_by_patch);
  void setup_internals();

private:
  std::vector<Particle> get_particles(int beg, int end);
  
public:
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

  const Grid_t& grid() const { return this->grid_; }

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

