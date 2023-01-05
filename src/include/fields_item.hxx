
#pragma once

#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "bnd.hxx"
#include "fields.hxx"
#include "fields3d.hxx"
#include "particles.hxx"
#include "psc_fields_c.h"

#include <mrc_profile.h>

#include <array>
#include <map>
#include <string>

// ======================================================================
// addKindSuffix

inline std::vector<std::string> addKindSuffix(
  const std::vector<std::string>& names, const Grid_t::Kinds& kinds)
{
  std::vector<std::string> result;
  for (int k = 0; k < kinds.size(); k++) {
    for (int m = 0; m < names.size(); m++) {
      result.emplace_back(names[m] + "_" + kinds[k].name);
    }
  }
  return result;
}

// ======================================================================
// ItemMomentBnd

template <typename S, typename Bnd>
class ItemMomentBnd
{
public:
  using storage_type = S;

  ItemMomentBnd(const Grid_t& grid) : bnd_{grid, grid.ibn} {}

  void add_ghosts(const Grid_t& grid, storage_type& mres_gt, const Int3& ib)
  {
    for (int p = 0; p < mres_gt.shape(4); p++) {
      add_ghosts_boundary(grid, mres_gt, ib, p, 0, mres_gt.shape(3));
    }

    bnd_.add_ghosts(grid, mres_gt, ib, 0, mres_gt.shape(3));
  }

private:
  // ----------------------------------------------------------------------
  // boundary stuff FIXME, should go elsewhere...

  template <typename FE>
  void add_ghosts_reflecting_lo(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                int p, int d, int mb, int me)
  {
    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iy = 0;
          {
            for (int m = mb; m < me; m++) {
              mres_gt(ix + ib[0], iy + ib[1], iz + ib[2], m, p) +=
                mres_gt(ix + ib[0], iy - 1 + ib[1], iz + ib[2], m, p);
            }
          }
        }
      }
    } else if (d == 2) {
      for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iz = 0;
          {
            for (int m = mb; m < me; m++) {
              mres_gt(ix + ib[0], iy + ib[1], iz + ib[2], m, p) +=
                mres_gt(ix + ib[0], iy + ib[1], iz - 1 + ib[2], m, p);
            }
          }
        }
      }
    } else {
      assert(0);
    }
  }

  template <typename FE>
  void add_ghosts_reflecting_hi(const Int3& ldims, FE& mres_gt, const Int3& ib,
                                int p, int d, int mb, int me)
  {
    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iy = ldims[1] - 1;
          {
            for (int m = mb; m < me; m++) {
              mres_gt(ix + ib[0], iy + ib[1], iz + ib[2], m, p) +=
                mres_gt(ix + ib[0], iy + 1 + ib[1], iz + ib[2], m, p);
            }
          }
        }
      }
    } else if (d == 2) {
      for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iz = ldims[2] - 1;
          {
            for (int m = mb; m < me; m++) {
              mres_gt(ix + ib[0], iy + ib[1], iz + ib[2], m, p) +=
                mres_gt(ix + ib[0], iy + ib[1], iz + 1 + ib[2], m, p);
            }
          }
        }
      }
    } else {
      assert(0);
    }
  }

  template <typename FE>
  void add_ghosts_boundary(const Grid_t& grid, FE& mres_gt, const Int3& ib,
                           int p, int mb, int me)
  {
    // lo
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryLo(p, d)) {
        if (grid.bc.prt_lo[d] == BND_PRT_REFLECTING ||
            grid.bc.prt_lo[d] == BND_PRT_OPEN) {
          add_ghosts_reflecting_lo(grid.ldims, mres_gt, ib, p, d, mb, me);
        }
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryHi(p, d)) {
        if (grid.bc.prt_hi[d] == BND_PRT_REFLECTING ||
            grid.bc.prt_hi[d] == BND_PRT_OPEN) {
          add_ghosts_reflecting_hi(grid.ldims, mres_gt, ib, p, d, mb, me);
        }
      }
    }
  }

private:
  Bnd bnd_;
};

// ======================================================================
// ItemMomentCRTP

template <typename Derived, typename S, typename Bnd = Bnd_>
class ItemMomentCRTP
{
public:
  using storage_type = S;
  using space_type = typename storage_type::space_type;
  using value_type = typename storage_type::value_type;

  static std::string name() { return Derived::moment_type::name(); }
  int n_comps() { return comp_names_.size(); }
  const std::vector<std::string>& comp_names() { return comp_names_; }

  explicit ItemMomentCRTP(const Grid_t& grid)
    : comp_names_{Derived::moment_type::comp_names(grid.kinds)}, bnd_{grid}
  {}

  template <typename Mparticles>
  auto operator()(const Mparticles& mprts)
  {
    Int3 ib = -mprts.grid().ibn;
    storage_type mres =
      psc::mflds::zeros<value_type, space_type>(mprts.grid(), n_comps(), ib);
    typename Derived::moment_type{}(mres, ib, mprts);
    bnd_.add_ghosts(mprts.grid(), mres, ib);
    return mres;
  }

protected:
  ItemMomentBnd<storage_type, Bnd> bnd_;

private:
  std::vector<std::string> comp_names_;
};
