
#pragma once

#include "injector_simple.hxx"
#include "const_accessor_simple.hxx"
#include "particles.hxx"
#include "particle_simple.hxx"
#include "particle_indexer.hxx"

#include <iterator>

// ======================================================================
// MparticlesPatchSimple

template<typename Mparticles>
struct MparticlesPatchSimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = typename Particle::real_t;
  using buf_t = std::vector<Particle>;
  using iterator = typename buf_t::iterator;
  using const_iterator = typename buf_t::const_iterator;

  // FIXME, I would like to delete the copy ctor because I don't
  // want to copy Patch by mistake, but that doesn't play well with
  // putting the patches into std::vector
  // MparticlesPatchSimple(const MparticlesPatchSimple&) = delete;

  MparticlesPatchSimple(Mparticles* mprts, int p)
    : mprts_(mprts),
      p_(p),
      grid_(&mprts->grid())
  {}

  Particle& operator[](int n) { return buf[n]; }
  const Particle& operator[](int n) const { return buf[n]; }
  const_iterator begin() const { return buf.begin(); }
  iterator begin() { return buf.begin(); }
  const_iterator end() const { return buf.end(); }
  iterator end() { return buf.end(); }
  unsigned int size() const { return buf.size(); }
  void reserve(unsigned int new_capacity) { buf.reserve(new_capacity); }

  void push_back(const Particle& new_prt)
  {
    // need to copy because we modify it
    auto prt = new_prt;
    checkInPatchMod(prt);
    validCellIndex(prt);
    buf.push_back(prt);
  }

  void resize(unsigned int new_size)
  {
    assert(new_size <= buf.capacity());
    buf.resize(new_size);
  }

  void check() const
  {
    for (auto& prt : buf) {
      validCellIndex(prt);
    }
  }

  // ParticleIndexer functionality
  int cellPosition(real_t xi, int d) const { return mprts_->pi_.cellPosition(xi, d); }
  int validCellIndex(const Particle& prt) const { return mprts_->pi_.validCellIndex(prt.x()); }

  void checkInPatchMod(Particle& prt) const { return mprts_->pi_.checkInPatchMod(prt.x()); }

  const Grid_t& grid() const { return *grid_; }
  const Mparticles& mprts() const { return *mprts_; }
  int p() const { return p_; }

  buf_t buf;

private:
  Mparticles* mprts_;
  int p_;
  const Grid_t* grid_;
};

// ======================================================================
// Mparticles

template<typename P>
struct MparticlesSimple : MparticlesBase
{
  using Particle = P;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;
  using Patch = MparticlesPatchSimple<MparticlesSimple>;
  using BndpParticle = P;
  using buf_t = typename Patch::buf_t;
  using Accessor = AccessorSimple<MparticlesSimple>;

  MparticlesSimple(const Grid_t& grid)
    : MparticlesBase(grid),
      pi_(grid)
  {
    patches_.reserve(grid.n_patches());
    for (int p = 0; p < grid.n_patches(); p++) {
      patches_.emplace_back(this, p);
    }
  }

  void reset(const Grid_t& grid) override
  {
    MparticlesBase::reset(grid);
    patches_.clear();
    patches_.reserve(grid.n_patches());
    for (int p = 0; p < grid.n_patches(); p++) {
      patches_.emplace_back(this, p);
    }
  }

  const Patch& operator[](int p) const { return patches_[p]; }
  Patch&       operator[](int p)       { return patches_[p]; }

  void reserve_all(const std::vector<uint> &n_prts_by_patch)
  {
    for (int p = 0; p < patches_.size(); p++) {
      patches_[p].reserve(n_prts_by_patch[p]);
    }
  }

  void resize_all(const std::vector<uint>& n_prts_by_patch)
  {
    for (int p = 0; p < patches_.size(); p++) {
      patches_[p].resize(n_prts_by_patch[p]);
    }
  }

  void reset() // FIXME, "reset" is used for two very different functions
  {
    for (int p = 0; p < patches_.size(); p++) {
      patches_[p].resize(0);
    }
  }

  std::vector<uint> get_size_all() const override
  {
    std::vector<uint> n_prts_by_patch(n_patches());
    for (int p = 0; p < patches_.size(); p++) {
      n_prts_by_patch[p] = patches_[p].size();
    }
    return n_prts_by_patch;
  }

  int get_n_prts() const override
  {
    int n_prts = 0;
    for (auto const& patch : patches_) {
      n_prts += patch.size();
    }
    return n_prts;
  }

  const ParticleIndexer<real_t>& particleIndexer() const { return pi_; }
  
  InjectorSimple<MparticlesSimple> injector() { return {*this}; }
  ConstAccessorSimple<MparticlesSimple> accessor() { return {*this}; }
  Accessor accessor_() { return {*this}; }
  
  void check() const
  {
    for (auto& patch: patches_) {
      patch.check();
    }
  }

  void dump(const std::string& filename)
  {
    FILE* file = fopen(filename.c_str(), "w");
    assert(file);

    for (int p = 0; p < n_patches(); p++) {
      auto& prts = (*this)[p];
      fprintf(file, "mparticles_dump: p%d n_prts = %d\n", p, prts.size());
      for (int n = 0; n < prts.size(); n++) {
	auto& prt = prts[n];
	fprintf(file, "mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g\n",
		n, prt.xi, prt.yi, prt.zi, prt.kind_,
		prt.pxi, prt.pyi, prt.pzi, prt.qni_wni_);
      }
    }
    fclose(file);
  }
  
  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  std::vector<Patch> patches_;
public: // FIXME
  ParticleIndexer<real_t> pi_;
};

