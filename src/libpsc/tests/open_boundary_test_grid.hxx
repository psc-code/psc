#pragma once

#include "grid.hxx"
#include "psc.h"

// ======================================================================
// makeOpenBoundaryTestGrid
//
// One patch of `gdims` unit cells, open (fields and particles) along the axes
// set in `open` and periodic along the others, with 2 ghosts along
// non-invariant axes. Electrons and ions have |q| = m = 1, and nicell = 1.

inline Grid_t* makeOpenBoundaryTestGrid(Int3 gdims, Vec3<bool> open, double dt)
{
  auto domain = Grid_t::Domain{gdims, Double3(gdims)};

  Int3 fld, prt, ibn;
  for (int d = 0; d < 3; d++) {
    fld[d] = open[d] ? BND_FLD_OPEN : BND_FLD_PERIODIC;
    prt[d] = open[d] ? BND_PRT_OPEN : BND_PRT_PERIODIC;
    ibn[d] = gdims[d] > 1 ? 2 : 0;
  }
  auto bc = psc::grid::BC{fld, fld, prt, prt};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;
  Grid_t::Normalization norm{norm_params};

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}
