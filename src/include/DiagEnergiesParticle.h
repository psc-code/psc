
#pragma once

#include <string>
#include <vector>

#include "psc_particles_double.h"

class DiagEnergiesParticle
{
public:
  std::vector<std::string> names() const { return {"E_electron", "E_ion"}; }

  template <typename Mparticles, typename MfieldsState>
  std::vector<double> operator()(Mparticles& mprts, MfieldsState& mflds) const
  {
    std::vector<double> vals(2);

    const auto& grid = mprts.grid();
    double fnqs = grid.norm.fnqs;
    double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];

    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      for (auto prt : accessor[p]) {
        double gamma =
          sqrt(1.f + sqr(prt.u()[0]) + sqr(prt.u()[1]) + sqr(prt.u()[2]));
        double Ekin = (gamma - 1.) * prt.m() * prt.w() * fnqs;
        double q = prt.q();
        if (q < 0.) {
          vals[0] += Ekin * fac;
        } else if (q > 0.) {
          vals[1] += Ekin * fac;
        } else {
          assert(0);
        }
      }
    }

    return vals;
  }
};
