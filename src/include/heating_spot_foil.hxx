
#pragma once

// ======================================================================
// HeatingSpotFoilParams

struct HeatingSpotFoilParams
{
  double zl; // in internal units (d_e)
  double zh;
  double xc;
  double yc;
  double rH;
  double T;
  double Mi;
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
    fac = (8.f * pow(T, 1.5)) / (sqrt(Mi) * width);
    // FIXME, I don't understand the sqrt(Mi) in here
  }
  
  double operator()(const double *crd)
  {
    double x = crd[0], y = crd[1], z = crd[2];

    if (z <= zl || z >= zh) {
      return 0;
    }
    
    return fac * (exp(-(sqr(x - (xc)) + sqr(y - (yc))) / sqr(rH)) +
		  exp(-(sqr(x - (xc)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
		  exp(-(sqr(x - (xc)) + sqr(y - (yc - Ly_))) / sqr(rH)) +
		  exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc))) / sqr(rH)) +
		  exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
		  exp(-(sqr(x - (xc + Lx_)) + sqr(y - (yc - Ly_))) / sqr(rH)) +
		  exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc))) / sqr(rH)) +
		  exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc + Ly_))) / sqr(rH)) +
		  exp(-(sqr(x - (xc - Lx_)) + sqr(y - (yc - Ly_))) / sqr(rH)));
  }

private:
  double fac;
  double Lx_, Ly_;
};

