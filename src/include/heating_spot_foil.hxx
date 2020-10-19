
#pragma once

// ======================================================================
// HeatingSpotFoilParams
#define HEATING_MAX_N_KINDS (10)
struct HeatingSpotFoilParams
{
  double zl; // in internal units (d_e)
  double zh;
  double xc;
  double yc;
  double rH;
  double T[HEATING_MAX_N_KINDS];
  double Mi;
  int n_kinds;
};

// ======================================================================
// HeatingSpotFoil

struct HeatingSpotFoil : HeatingSpotFoilParams
{
  HeatingSpotFoil() = default;

  HeatingSpotFoil(const Grid_t& grid, const HeatingSpotFoilParams& params)
    : HeatingSpotFoilParams(params),
      Lx_(grid.domain.length[0]),
      Ly_(grid.domain.length[1])
  {
    double width = zh - zl;

    assert(n_kinds < HEATING_MAX_N_KINDS);
    // inialize a fac for each population
    for (int i = 0; i < n_kinds; i++)
      fac[i] = (8.f * pow(T[i], 1.5)) / (sqrt(Mi) * width);

    // FIXME, I don't understand the sqrt(Mi) in here
  }

  template <typename R>
  KG_INLINE R operator()(const R* crd, int kind)
  {
    R x = crd[0], y = crd[1], z = crd[2];
    if (fac[kind] == 0.0)
      return 0;
    if (z <= zl || z >= zh) {
      return 0;
    }

    return fac[kind] *
           (std::exp(-(sqr(x - (xc)) + sqr(y - (yc))) / sqr(rH))); // +
#if 0
            std::exp(-(sqr(x - (xc)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc)) + sqr(y - (yc - Ly_))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc - Ly_))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
            std::exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc - Ly_))) / sqr(rH)));
#endif
  }

  // private:
  double fac[HEATING_MAX_N_KINDS];
  double Lx_, Ly_;
};
