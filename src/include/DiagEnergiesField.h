
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

  std::vector<double> operator()(MparticlesBase& mprts,
                                 MfieldsStateBase& mflds_base) const
  {
    std::vector<double> EH2(6);

    const Grid_t& grid = mprts.grid();
    auto& mf = mflds_base.get_as<MfieldsStateDouble>(EX, HX + 3);
    for (int p = 0; p < grid.n_patches(); p++) {
      double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];
      auto F = mf[p];
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
    mflds_base.put_as(mf, 0, 0);
    return EH2;
  }
};
