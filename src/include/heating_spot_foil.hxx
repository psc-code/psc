
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
  
  HeatingSpotFoil(const HeatingSpotFoilParams& params)
    : HeatingSpotFoilParams(params)
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
    
    return fac * exp(-(sqr(x-xc) + sqr(y-yc)) / sqr(rH));
  }

private:
  double fac;
};

