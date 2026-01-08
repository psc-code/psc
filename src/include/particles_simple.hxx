
#pragma once

#include "injector_simple.hxx"
#include "const_accessor_simple.hxx"
#include "particles.hxx"
#include "particle_simple.hxx"
#include "particle_indexer.hxx"
#include "UniqueIdGenerator.h"

#include <iterator>

// ======================================================================

template <typename T>
class Span
{
public:
  using element_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator = T*;
  using const_iterator = const T*;

  Span(pointer data, size_t size) : begin_{data}, end_{data + size} {}

  iterator begin() const { return begin_; }
  iterator end() const { return end_; }
  size_t size() const { return end_ - begin_; }

private:
  pointer begin_;
  pointer end_;
};

// ======================================================================
// MparticlesStorage

template <typename _Particle>
struct MparticlesStorage
{
  using Particle = _Particle;
  using PatchBuffer = std::vector<Particle>;
  using Buffers = std::vector<PatchBuffer>;
  using Range = Span<Particle>;
  using iterator = typename Range::iterator;
  using const_iterator = typename Range::const_iterator;

  MparticlesStorage(uint n_patches) : bufs_(n_patches) {}

  void reset(const Grid_t& grid) { bufs_ = Buffers(grid.n_patches()); }

  void reserve_all(const std::vector<uint>& n_prts_by_patch)
  {
    for (int p = 0; p < bufs_.size(); p++) {
      bufs_[p].reserve(n_prts_by_patch[p]);
    }
  }

  void resize_all(const std::vector<uint>& n_prts_by_patch)
  {
    for (int p = 0; p < bufs_.size(); p++) {
      assert(n_prts_by_patch[p] <= bufs_[p].capacity());
      bufs_[p].resize(n_prts_by_patch[p]);
    }
  }

  void clear()
  {
    for (int p = 0; p < bufs_.size(); p++) {
      bufs_[p].resize(0);
    }
  }

  std::vector<uint> sizeByPatch() const
  {
    std::vector<uint> n_prts_by_patch(bufs_.size());
    for (int p = 0; p < bufs_.size(); p++) {
      n_prts_by_patch[p] = bufs_[p].size();
    }
    return n_prts_by_patch;
  }

  int size() const
  {
    int n_prts = 0;
    for (const auto& buf : bufs_) {
      n_prts += buf.size();
    }
    return n_prts;
  }

  Range operator[](int p) { return {bufs_[p].data(), bufs_[p].size()}; }
  Particle& at(int p, int n)
  {
    return bufs_[p][n];
  } // FIXME, ugly and not great for effciency
  void push_back(int p, const Particle& prt) { bufs_[p].push_back(prt); }

  Buffers& bndBuffers() { return bufs_; }

private:
  Buffers bufs_;
};

// ======================================================================
// MparticlesSimple

template <typename P>
struct MparticlesSimple : MparticlesBase
{
  using Particle = P;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;
  using BndpParticle = P;
  using Accessor = AccessorSimple<MparticlesSimple>;
  using ConstAccessor = ConstAccessorSimple<MparticlesSimple>;
  using Storage = MparticlesStorage<Particle>;
  using BndBuffer = typename Storage::PatchBuffer;
  using BndBuffers = typename Storage::Buffers;

  struct Patch
  {
    using iterator = typename Storage::iterator;
    using const_iterator = typename Storage::const_iterator;

    Patch(MparticlesSimple& mprts, int p) : mprts_(mprts), p_(p) {}

    Patch(const Patch&) = delete;
    Patch(Patch&&) = default;

    Particle& operator[](int n) { return mprts_.storage_.at(p_, n); }
    const Particle& operator[](int n) const
    {
      return mprts_.storage_.at(p_, n);
    }

    iterator begin() { return mprts_.storage_[p_].begin(); }
    iterator end() { return mprts_.storage_[p_].end(); }
    unsigned int size() const { return mprts_.storage_[p_].size(); }

    const Grid_t& grid() const { return mprts_.grid(); }
    const MparticlesSimple& mprts() const { return mprts_; }
    int p() const { return p_; }

  private:
    MparticlesSimple& mprts_;
    int p_;
  };

  explicit MparticlesSimple(const Grid_t& grid)
    : MparticlesBase(grid),
      storage_(grid.n_patches()),
      uid_gen(grid.comm()),
      pi_(grid)
  {}

  MparticlesSimple(const MparticlesSimple&) = delete;
  MparticlesSimple(MparticlesSimple&& o) = default;

  MparticlesSimple& operator=(MparticlesSimple&& o) = default;

  void reset(const Grid_t& grid) override
  {
    MparticlesBase::reset(grid);
    storage_.reset(grid);
  }

  Patch operator[](int p) const
  {
    return {const_cast<MparticlesSimple&>(*this), p};
  } // FIXME, isn't actually const

  void reserve_all(const std::vector<uint>& n_prts_by_patch)
  {
    storage_.reserve_all(n_prts_by_patch);
  }
  void resize_all(const std::vector<uint>& n_prts_by_patch)
  {
    storage_.resize_all(n_prts_by_patch);
  }
  void clear() { storage_.clear(); }
  std::vector<uint> sizeByPatch() const override
  {
    return storage_.sizeByPatch();
  }
  int size() const override { return storage_.size(); }

  const ParticleIndexer<real_t>& particleIndexer() const { return pi_; }

  int cellPosition(int p, int n, int d) const
  {
    return pi_.cellPosition((*this)[p][n].x[d], d);
  }

  int validCellIndex(const Real3& x) const { return pi_.validCellIndex(x); }

  void checkInPatchMod(Particle& prt) const
  {
    return pi_.checkInPatchMod(prt.x);
  }

  void push_back(int p, const Particle& new_prt)
  {
    // need to copy because we modify it
    auto prt = new_prt;
    checkInPatchMod(prt);
    validCellIndex(prt.x);
    storage_.push_back(p, prt);
  }

  InjectorSimple<MparticlesSimple> injector() { return {*this}; }
  ConstAccessor accessor() const
  {
    return {const_cast<MparticlesSimple&>(*this)};
  } // FIXME
  Accessor accessor_() { return {*this}; }

  BndBuffers& bndBuffers() { return storage_.bndBuffers(); }

  void dump(const std::string& filename)
  {
    FILE* file = fopen(filename.c_str(), "w");
    assert(file);

    for (int p = 0; p < n_patches(); p++) {
      auto& prts = (*this)[p];
      fprintf(file, "mparticles_dump: p%d n_prts = %d\n", p, prts.size());
      for (int n = 0; n < prts.size(); n++) {
        auto& prt = prts[n];
        fprintf(file,
                "mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g\n", n,
                prt.xi, prt.yi, prt.zi, prt.kind_, prt.pxi, prt.pyi, prt.pzi,
                prt.qni_wni_);
      }
    }
    fclose(file);
  }

  real_t prt_q(const Particle& prt) const { return grid().kinds[prt.kind].q; }

  real_t prt_m(const Particle& prt) const { return grid().kinds[prt.kind].m; }

  real_t prt_w(const Particle& prt) const { return prt.qni_wni / prt_q(prt); }

  void define_species(const char* name, double q, double m, double max_local_np,
                      double max_local_nm, double sort_interval,
                      double sort_out_of_place)
  {}

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  Storage storage_;

public: // FIXME
  psc::particle::UniqueIdGenerator uid_gen;
  ParticleIndexer<real_t> pi_;
};

template <>
const MparticlesSimple<ParticleSimple<float>>::Convert
  MparticlesSimple<ParticleSimple<float>>::convert_to_;
extern template const MparticlesSimple<ParticleSimple<float>>::Convert
  MparticlesSimple<ParticleSimple<float>>::convert_to_;
template <>
const MparticlesSimple<ParticleSimple<float>>::Convert
  MparticlesSimple<ParticleSimple<float>>::convert_from_;
extern template const MparticlesSimple<ParticleSimple<float>>::Convert
  MparticlesSimple<ParticleSimple<float>>::convert_from_;

template <>
const MparticlesSimple<ParticleSimple<double>>::Convert
  MparticlesSimple<ParticleSimple<double>>::convert_to_;
extern template const MparticlesSimple<ParticleSimple<double>>::Convert
  MparticlesSimple<ParticleSimple<double>>::convert_to_;
template <>
const MparticlesSimple<ParticleSimple<double>>::Convert
  MparticlesSimple<ParticleSimple<double>>::convert_from_;
extern template const MparticlesSimple<ParticleSimple<double>>::Convert
  MparticlesSimple<ParticleSimple<double>>::convert_from_;
