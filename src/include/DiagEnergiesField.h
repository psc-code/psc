
#pragma once

#include <string>
#include <vector>

#include "psc_fields_c.h"
#include "fields.hxx"

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
    auto&& h_gt_mflds = gt::host_mirror(mflds.storage());
    gt::copy(mflds.storage(), h_gt_mflds);

    double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];
    for (int p = 0; p < grid.n_patches(); p++) {
      auto bnd = mflds.ibn();
      auto flds = h_gt_mflds.view(_s(-bnd[0], bnd[0]), _s(-bnd[1], bnd[1]),
                                  _s(-bnd[2], bnd[2]), _all, p);
      // FIXME, this doesn't handle non-periodic b.c. right
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        EH2[0] += sqr(flds(i, j, k, EX)) * fac;
        EH2[1] += sqr(flds(i, j, k, EY)) * fac;
        EH2[2] += sqr(flds(i, j, k, EZ)) * fac;
        EH2[3] += sqr(flds(i, j, k, HX)) * fac;
        EH2[4] += sqr(flds(i, j, k, HY)) * fac;
        EH2[5] += sqr(flds(i, j, k, HZ)) * fac;
      });
    }
    return EH2;
  }
};
