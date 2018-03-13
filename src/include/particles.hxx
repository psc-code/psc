
#ifndef PARTICLES_HXX
#define PARTICLES_HXX

#include "psc_particles.h"
#include "grid.hxx"
#include "psc_bits.h"

#include <mrc_profile.h>

#include <vector>
#include <unordered_map>
#include <typeindex>
#include <string>

// ======================================================================
// MparticlesBase

struct particle_base_t
{
  struct real_t {};
};

struct patch_base_t
{
  struct buf_t {};
};

struct MparticlesBase
{
  using particle_t = particle_base_t;
  using patch_t = patch_base_t;

  using copy_func_t = void (*)(MparticlesBase&, MparticlesBase&);
  using Convert = std::unordered_map<std::type_index, copy_func_t>;
  
  MparticlesBase(const Grid_t& grid)
    : grid_(&grid)
  {}

  const Grid_t& grid() const { return *grid_; }
  int n_patches() const { return grid_->n_patches(); }

  virtual ~MparticlesBase() {}
  virtual int get_n_prts() const = 0;
  virtual void get_size_all(uint *n_prts_by_patch) const = 0;
  virtual void reserve_all(const uint *n_prts_by_patch) = 0;
  virtual void resize_all(const uint *n_prts_by_patch) = 0;
  virtual void inject(int p, const psc_particle_inject& new_prt) { assert(0); }
  virtual void inject_reweight(int p, const psc_particle_inject& new_prt) { assert(0); }
  virtual MparticlesBase* create(const Grid_t& grid) = 0;

  virtual const Convert& convert_to() { static const Convert convert_to_; return convert_to_; }
  virtual const Convert& convert_from() { static const Convert convert_from_; return convert_from_; }
  static void convert(MparticlesBase& mp_from, MparticlesBase& mp_to);
  
protected:
  const Grid_t* grid_;
public:
  bool inited = true; // FIXME hack to avoid dtor call when not yet constructed
};

using psc_mparticles_copy_func_t = MparticlesBase::copy_func_t; // FIXME, get rid of

// ======================================================================
// PscMparticles

template<typename S>
struct PscMparticles;

using PscMparticlesBase = PscMparticles<MparticlesBase>;

PscMparticlesBase PscMparticlesCreate(MPI_Comm comm, const Grid_t& grid, const char *type);

template<typename S>
struct PscMparticles
{
  using Self = PscMparticles<S>;
  using sub_t = S;
  using particle_t = typename sub_t::particle_t;
  using real_t = typename particle_t::real_t;
  using patch_t = typename sub_t::patch_t;
  using particle_buf_t = typename patch_t::buf_t;

  explicit PscMparticles(psc_mparticles *mprts)
    : mprts_(mprts)
  {
    if (mprts != nullptr) {
      assert(dynamic_cast<sub_t*>(mrc_to_subobj(mprts, MparticlesBase)));
    }
  }

  static Self create(MPI_Comm comm, const Grid_t& grid)
  {
    auto mprts = PscMparticlesCreate(comm, grid, mparticles_traits<Self>::name);
    return Self{mprts.mprts()}; // odd way of returning derived type
  }
  
  template<typename MP>
  MP get_as(uint flags = 0);

  template<typename MP>
  void put_as(MP mprts_base, uint flags = 0);
  
  psc_mparticles *mprts() { return mprts_; }
  
  sub_t* operator->() { return sub(); }

  sub_t* sub() { return mrc_to_subobj(mprts_, sub_t); }

  patch_t& operator[](int p)
  {
    return (*this->sub())[p];
  }

private:
  static void convert(MparticlesBase& mp_from, MparticlesBase& mp_to);

private:
  psc_mparticles *mprts_;
};

inline PscMparticlesBase PscMparticlesCreate(MPI_Comm comm, const Grid_t& grid, const char *type)
{
  psc_mparticles* mprts = psc_mparticles_create(comm);
  psc_mparticles_set_type(mprts, type);
  mprts->grid = &grid;
  psc_mparticles_setup(mprts);
  return PscMparticlesBase{mprts};
}

// ======================================================================
// MparticlesWrapper

template<typename Mparticles>
class MparticlesWrapper
{
public:
  const static size_t size = sizeof(Mparticles);

  constexpr static const char* name = mparticles_traits<PscMparticles<Mparticles>>::name;
  
  static void setup(struct psc_mparticles* _mprts)
  {
    new(_mprts->obj.subctx) Mparticles{*_mprts->grid};
  }

  static void destroy(struct psc_mparticles* _mprts)
  {
    if (!mrc_to_subobj(_mprts, MparticlesBase)->inited) return; // FIXME
    PscMparticles<Mparticles> mprts(_mprts);
    mprts->~Mparticles();
  }
};

template<typename Mparticles>
struct psc_mparticles_ops_ : psc_mparticles_ops {
  using Wrapper_t = MparticlesWrapper<Mparticles>;
  psc_mparticles_ops_() {
    name    = Wrapper_t::name;
    size    = Wrapper_t::size;
    setup   = Wrapper_t::setup;
    destroy = Wrapper_t::destroy;
  }

  MparticlesBase* create(const Grid_t&grid) { return new Mparticles{grid}; }
};

#endif

