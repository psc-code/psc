
#pragma once

#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "bnd.hxx"
#include "centering.hxx"
#include "fields.hxx"
#include "fields3d.hxx"
#include "particles.hxx"
#include "psc_fields_c.h"
#include "add_ghosts_reflecting.hxx"

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

template <typename FE>
void add_ghosts_boundary(const Grid_t& grid, FE& mres_gt, const Int3& ib, int p,
                         int mb, int me, centering::Centering c)
{
  if (c == centering::CC) {
    // lo
    for (int d = 0; d < 3; d++) {
      // FIXME why reflect for open BCs?
      if (grid.atBoundaryLo(p, d) && (grid.bc.prt_lo[d] == BND_PRT_REFLECTING ||
                                      grid.bc.prt_lo[d] == BND_PRT_OPEN)) {
        add_ghosts_reflecting_lo_cc(grid.ldims, mres_gt, ib, p, d, mb, me);
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      // FIXME why reflect for open BCs?
      if (grid.atBoundaryHi(p, d) && (grid.bc.prt_hi[d] == BND_PRT_REFLECTING ||
                                      grid.bc.prt_hi[d] == BND_PRT_OPEN)) {
        add_ghosts_reflecting_hi_cc(grid.ldims, mres_gt, ib, p, d, mb, me);
      }
    }
  } else if (c == centering::NC) {
    // lo
    for (int d = 0; d < 3; d++) {
      // FIXME why reflect for open BCs?
      if (grid.atBoundaryLo(p, d) && (grid.bc.prt_lo[d] == BND_PRT_REFLECTING ||
                                      grid.bc.prt_lo[d] == BND_PRT_OPEN)) {
        add_ghosts_reflecting_lo_nc(grid.ldims, mres_gt, ib, p, d, mb, me);
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      // FIXME why reflect for open BCs?
      if (grid.atBoundaryHi(p, d) && (grid.bc.prt_hi[d] == BND_PRT_REFLECTING ||
                                      grid.bc.prt_hi[d] == BND_PRT_OPEN)) {
        add_ghosts_reflecting_hi_nc(grid.ldims, mres_gt, ib, p, d, mb, me);
      }
    }
  }
}

// ======================================================================
// ItemMomentBnd

template <typename S, typename Bnd, centering::Centering C>
class ItemMomentBnd
{
public:
  using storage_type = S;

  ItemMomentBnd(const Grid_t& grid) {}

  void add_ghosts(const Grid_t& grid, storage_type& mres_gt, const Int3& ib)
  {
    for (int p = 0; p < mres_gt.shape(4); p++) {
      add_ghosts_boundary(grid, mres_gt, ib, p, 0, mres_gt.shape(3), C);
    }

    bnd_.add_ghosts(grid, mres_gt, ib, 0, mres_gt.shape(3));
  }

private:
  Bnd bnd_;
};

// ======================================================================
// ItemMoment
//
// This build upon psc::mmoment::moment_* and
// * stores comp_names, provides n_comps()
// * allocates the result of the moment calculation
// * performs add_ghosts

template <typename MT, typename S, typename Bnd = Bnd_>
class ItemMoment
{
public:
  using moment_type = MT;
  using storage_type = S;
  using space_type = typename storage_type::space_type;
  using value_type = typename storage_type::value_type;

  static std::string name() { return moment_type::name(); }
  int n_comps() { return comp_names_.size(); }
  const std::vector<std::string>& comp_names() { return comp_names_; }

  explicit ItemMoment(const Grid_t& grid)
    : comp_names_{moment_type::comp_names(grid.kinds)}, bnd_{grid}
  {}

  template <typename Mparticles>
  auto operator()(const Mparticles& mprts)
  {
    Int3 ib = -mprts.grid().ibn;
    // FIXME, gt::gtensor and psc::gtensor are slightly different, and ideally
    // zeros() shouldn't actually allocate, but probably it does, so this
    // wastes memory and a copy
    storage_type mres =
      psc::mflds::zeros<value_type, space_type>(mprts.grid(), n_comps(), ib);
    moment_type{}(mres, ib, mprts);
    bnd_.add_ghosts(mprts.grid(), mres, ib);
    return mres;
  }

private:
  ItemMomentBnd<storage_type, Bnd, moment_type::CENTERING> bnd_;
  std::vector<std::string> comp_names_;
};
