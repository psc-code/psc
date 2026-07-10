#include "kg/Vec3.h"

template <typename real_t>
struct RadiatingBoundary
{
  using Real3 = Vec3<real_t>;

  virtual real_t sample_exterior_field(int m, double t, int p, Real3 x3) = 0;

  virtual void tick(double t) {}
};