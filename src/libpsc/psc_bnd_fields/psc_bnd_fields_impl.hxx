#pragma once

#include "psc.h"
#include "kg/VecRange.hxx"
#include "fields.hxx"
#include "bnd_fields.hxx"
#include "radiating_bnd.hxx"
#include "field_bc_util.hxx"

#include <mrc_bits.h>

#include <limits>

// #define DEBUG

template <typename MFIELDS_STATE, typename Dim>
struct BndFields_ : BndFieldsBase
{
  using Self = BndFields_<MFIELDS_STATE, Dim>;
  using MfieldsState = MFIELDS_STATE;
  using real_t = typename MfieldsState::real_t;
  using Real3 = Vec3<real_t>;
  using fields_view_t = typename MfieldsState::fields_view_t;
  using dim_t = Dim;

  // ----------------------------------------------------------------------
  // fill_ghosts_E

  void fill_ghosts_E(MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();

    for (int p = 0; p < mflds.n_patches(); p++) {
      // lo
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryLo(p, d)) {
          switch (grid.bc.fld_lo[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }

      // hi
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryHi(p, d)) {
          switch (grid.bc.fld_hi[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // fill_ghosts_H

  void fill_ghosts_H(MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();

    if (radiation) {
    }

    for (int p = 0; p < mflds.n_patches(); p++) {
      // lo
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryLo(p, d)) {
          switch (grid.bc.fld_lo[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
      // hi
      for (int d = 0; d < 3; d++) {
        if (grid.atBoundaryHi(p, d)) {
          switch (grid.bc.fld_hi[d]) {
            case BND_FLD_PERIODIC: {
              break;
            }
            case BND_FLD_CONDUCTING_WALL: {
              break;
            }
            case BND_FLD_OPEN: {
              break;
            }
            default: {
              assert(0);
            }
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // add_ghosts_J

  void add_ghosts_J(MfieldsState& mflds) {}

  RadiatingBoundary<real_t>* radiation = nullptr;
};

// ======================================================================
// BndFieldsNone

// used by CUDA

template <class MFIELDS_STATE>
struct BndFieldsNone : BndFieldsBase
{
  using MfieldsState = MFIELDS_STATE;

  // clang-format off
  void fill_ghosts_E(MfieldsState& mflds) {};
  void fill_ghosts_H(MfieldsState& mflds) {};
  void add_ghosts_J(MfieldsState& mflds) {};
  // clang-format on
};
