
#pragma once

#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles_sort.hxx"
#include "injector_buffered.hxx"
#include "mparticles_patch_cuda.hxx"

#include "particles.hxx"
#include "psc_bits.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/iterator/zip_iterator.h>

#include <gtensor/span.h>

extern std::size_t mem_particles;

template <typename T>
struct MparticlesCudaStorage_;

using MparticlesCudaStorage =
  MparticlesCudaStorage_<psc::device_vector<float4>>;

using HMparticlesCudaStorage =
  MparticlesCudaStorage_<thrust::host_vector<float4>>;

using DMparticlesCudaStorage = MparticlesCudaStorage_<gt::span<float4>>;

static const double growth_factor = 1.5;

// ======================================================================
// helper

namespace helper
{
template <typename S>
void reserve(S& vec, std::size_t n)
{
  assert(n >= vec.size());
  S new_vec(n);
  new_vec.resize(vec.size());
  thrust::copy(vec.begin(), vec.end(), new_vec.begin());
  vec = std::move(new_vec);
}
} // namespace helper

// ======================================================================
// MparticlesCudaStorage_

template <typename S>
class MparticlesCudaStorage_
{
  using xi4_iterator = typename S::iterator;
  using pxi4_iterator = typename S::iterator;
  using iterator_tuple = thrust::tuple<xi4_iterator, pxi4_iterator>;

public:
  using value_type = typename S::value_type[2];
  using iterator = thrust::zip_iterator<iterator_tuple>;

  MparticlesCudaStorage_() = default;
  MparticlesCudaStorage_(const MparticlesCudaStorage_& other) = default;
  MparticlesCudaStorage_& operator=(const MparticlesCudaStorage_& other) =
    delete;

  MparticlesCudaStorage_(uint n) : xi4(n), pxi4(n)
  {
    mem_add(xi4);
    mem_add(pxi4);
  }

  MparticlesCudaStorage_(const S& xi4, const S& pxi4) : xi4{xi4}, pxi4{pxi4} {}

  template <typename U>
  GT_INLINE void mem_add(const U&)
  {}

  template <typename T>
  GT_INLINE void mem_add(const psc::device_vector<T>& v)
  {
#ifndef __CUDA_ARCH__
    mem_particles += allocated_bytes(v);
#endif
  }

  template <typename U>
  GT_INLINE void mem_sub(const U&)
  {}

  template <typename T>
  GT_INLINE void mem_sub(const psc::device_vector<T>& v)
  {
#ifndef __CUDA_ARCH__
    mem_particles -= allocated_bytes(v);
#endif
  }

  GT_INLINE ~MparticlesCudaStorage_()
  {
    mem_sub(xi4);
    mem_sub(pxi4);
  }

  template <typename SO>
  MparticlesCudaStorage_& operator=(const MparticlesCudaStorage_<SO>& other)
  {
    mem_sub(xi4);
    xi4 = other.xi4;
    mem_add(xi4);

    mem_sub(pxi4);
    pxi4 = other.pxi4;
    mem_add(pxi4);
    return *this;
  }

  template <typename SO>
  __host__ MparticlesCudaStorage_(const SO& other)
    : xi4{other.xi4}, pxi4{other.pxi4}
  {
    mem_add(xi4);
    mem_add(pxi4);
  }

  template <typename IT>
  MparticlesCudaStorage_(IT first, IT last)
  {
    resize(last - first);
    thrust::copy(first, last, begin());
  }

  __host__ __device__ size_t size() const { return xi4.size(); }

  __host__ __device__ iterator begin()
  {
    return iterator({xi4.begin(), pxi4.begin()});
  }

  __host__ __device__ iterator end()
  {
    return iterator({xi4.end(), pxi4.end()});
  }

  __host__ void resize(size_t n)
  {
    mem_sub(xi4);
    mem_sub(pxi4);

    // grow arrays by fixed growth factor only
    if (n > xi4.capacity()) {
#if 0
      xi4.reserve(growth_factor * n);
      pxi4.reserve(growth_factor * n);
#else
      // work around thrust mem leak and exponential growth
      helper::reserve(xi4, growth_factor * n);
      helper::reserve(pxi4, growth_factor * n);
#endif
    }
    xi4.resize(n);
    pxi4.resize(n);

    mem_add(xi4);
    mem_add(pxi4);
  }

  __host__ __device__ DParticleCuda operator[](int n) const
  {
    float4 xi4 = this->xi4[n], pxi4 = this->pxi4[n];
    return DParticleCuda{
      {xi4.x, xi4.y, xi4.z},    {pxi4.x, pxi4.y, pxi4.z}, pxi4.w,
      cuda_float_as_int(xi4.w), psc::particle::Id(),      psc::particle::Tag()};
  }

  __host__ __device__ void store(const DParticleCuda& prt, int n)
  {
    store_position(prt, n);
    store_momentum(prt, n);
  }

  __host__ __device__ void store_position(const DParticleCuda& prt, int n)
  {
    xi4[n] = {prt.x[0], prt.x[1], prt.x[2], cuda_int_as_float(prt.kind)};
  }

  __host__ __device__ void store_momentum(const DParticleCuda& prt, int n)
  {
    pxi4[n] = {prt.u[0], prt.u[1], prt.u[2], prt.qni_wni};
  }

  friend void swap(MparticlesCudaStorage_& first,
                   MparticlesCudaStorage_& second)
  {
    using std::swap;
    swap(first.xi4, second.xi4);
    swap(first.pxi4, second.pxi4);
  }

  DMparticlesCudaStorage to_kernel()
  {
    return DMparticlesCudaStorage{{xi4.data().get(), xi4.size()},
                                  {pxi4.data().get(), pxi4.size()}};
  }

  S xi4;
  S pxi4;
};

// ======================================================================
// cuda_mparticles_base

template <typename BS, typename S>
struct cuda_mparticles_base : cuda_mparticles_indexer<BS>
{
  using storage_type = S;

  cuda_mparticles_base(const Grid_t& grid);
  // copy constructor would work fine, but we don't want to copy everything
  // by accident
  cuda_mparticles_base(const cuda_mparticles<BS>&) = delete;

protected:
  void resize(uint size);

public:
  std::vector<uint> sizeByPatch() const;
  void clear()
  {
    storage.resize(0);
    by_block_.clear();
    n_prts = 0;
  }

  // per particle
  storage_type storage;

  // per block
  cuda_mparticles_sort_by_block by_block_;

  uint n_prts = 0; // total # of particles across all patches
  const Grid_t& grid_;
};

// ======================================================================
// cuda_mparticles

template <typename BS>
struct DMparticlesCuda;

template <typename _BS>
struct cuda_mparticles : cuda_mparticles_base<_BS, MparticlesCudaStorage>
{
  using base_type = cuda_mparticles_base<_BS, MparticlesCudaStorage>;
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
  ~cuda_mparticles();

  uint size();
  void inject(const std::vector<Particle>& buf,
              const std::vector<uint>& buf_n_by_patch);

  std::vector<uint> get_offsets() const;
  std::vector<Particle> get_particles();

  double mem_fraction() const
  {
    return double(this->storage.xi4.size()) / this->storage.xi4.capacity();
  }

  void dump(const std::string& filename) const;
  void dump_by_patch(uint* n_prts_by_patch);

  // internal / testing use
  void inject_initial(const std::vector<Particle>& buf,
                      const std::vector<uint>& n_prts_by_patch);
  void setup_internals();

public:
  void reorder();
  void reorder(const psc::device_vector<uint>& d_id);
  void reorder_and_offsets(const psc::device_vector<uint>& d_idx,
                           const psc::device_vector<uint>& d_id,
                           psc::device_vector<uint>& d_off);
  void reorder_and_offsets_slow();
  void swap_alt();

  bool check_in_patch_unordered_slow();
  bool check_bidx_id_unordered_slow();
  bool check_ordered();
  bool check_bidx_after_push();

  void resize(uint n_prts);

  const Grid_t& grid() const { return this->grid_; }

public:
  MparticlesCudaStorage
    alt_storage; // storage for out-of-place reordering of particle data

  std::vector<Real3> xb_by_patch; // lower left corner for each patch

  bool need_reorder = {
    false}; // particles haven't yet been put into their sorted order
};

template <typename BS_>
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
      storage{cmprts.storage.to_kernel()},
      alt_storage{cmprts.alt_storage.to_kernel()},
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

  __device__ real_t prt_q(const DParticleCuda& prt) const
  {
    return q(prt.kind);
  }

  __device__ real_t prt_m(const DParticleCuda& prt) const
  {
    return m(prt.kind);
  }

  __device__ real_t prt_w(const DParticleCuda& prt) const
  {
    return prt.qni_wni / prt_q(prt);
  }

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
  uint* off_;
  uint* bidx_;
  uint* id_;
  uint n_blocks_;
};
