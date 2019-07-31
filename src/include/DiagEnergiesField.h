
#pragma once

#include <string>
#include <vector>

#include "psc_fields_c.h"

class DiagEnergiesField
{
public:
  std::vector<std::string> names() const
  {
    return {"EX2", "EY2", "EZ2", "BX2", "BY2", "BZ2"};
  }

  template <typename Mparticles, typename MfieldsState>
  std::vector<double> operator()(Mparticles& mprts, MfieldsState& mflds) const
  {
    std::vector<double> EH2(6);

    const Grid_t& grid = mprts.grid();
    for (int p = 0; p < grid.n_patches(); p++) {
      double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];
      auto F = mflds[p];
      // FIXME, this doesn't handle non-periodic b.c. right
      grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
        EH2[0] += sqr(F(EX, ix, iy, iz)) * fac;
        EH2[1] += sqr(F(EY, ix, iy, iz)) * fac;
        EH2[2] += sqr(F(EZ, ix, iy, iz)) * fac;
        EH2[3] += sqr(F(HX, ix, iy, iz)) * fac;
        EH2[4] += sqr(F(HY, ix, iy, iz)) * fac;
        EH2[5] += sqr(F(HZ, ix, iy, iz)) * fac;
      });
    }

    return EH2;
  }
};
