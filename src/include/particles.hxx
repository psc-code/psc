
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

struct MparticlesBase
{
  using convert_func_t = void (*)(MparticlesBase&, MparticlesBase&);
  using Convert = std::unordered_map<std::type_index, convert_func_t>;

  MparticlesBase(const Grid_t& grid) : grid_(&grid) {}

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }

  const Grid_t& grid() const { return *grid_; }
  int n_patches() const { return grid_->n_patches(); }

  virtual ~MparticlesBase() {}
  virtual int size() const = 0;
  virtual std::vector<uint> sizeByPatch() const = 0;

  template <typename MP>
  MP& get_as(uint flags = 0)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(*this) == typeid(MP)) {
      return *dynamic_cast<MP*>(this);
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mparticles_get_as", 1., 0, 0);
    }

    prof_start(pr);
    auto& mprts = *new MP{grid()};
    assert(!(flags & MP_DONT_COPY));
    MparticlesBase::convert(*this, mprts);
    prof_stop(pr);
    return mprts;
  }

  template <typename MP>
  void put_as(MP& mprts, uint flags = 0)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(mprts) == typeid(*this)) {
      return;
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mparticles_put_as", 1., 0, 0);
    }
    prof_start(pr);
    // mprintf("put_as %s -> %s from\n", typeid(mprts).name(),
    // typeid(*thi).name());

    if (!(flags & MP_DONT_COPY)) {
      MparticlesBase::convert(mprts, *this);
    }

    // mprintf("put_as %s -> %s to\n", typeid(mprts).name(),
    // typeid(*this).name());
    delete &mprts;

    prof_stop(pr);
  }

  void view() const
  {
    MPI_Comm comm = MPI_COMM_WORLD; // FIXME
    mpi_printf(comm, "  n_patches    = %d\n", n_patches());
    mpi_printf(comm, "  n_prts_total = %d\n", size());

    auto n_prts_by_patch = sizeByPatch();

    for (int p = 0; p < n_patches(); p++) {
      mpi_printf(comm, "  p %d: n_prts = %d\n", p, n_prts_by_patch[p]);
    }
  }

  virtual const Convert& convert_to()
  {
    static const Convert convert_to_;
    return convert_to_;
  }
  virtual const Convert& convert_from()
  {
    static const Convert convert_from_;
    return convert_from_;
  }
  static void convert(MparticlesBase& mp_from, MparticlesBase& mp_to);

protected:
  const Grid_t* grid_;
};

template <typename Particle>
class ForComponents;

#endif
