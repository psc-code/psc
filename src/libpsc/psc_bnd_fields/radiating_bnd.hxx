#include "kg/Vec3.h"

template <typename real_t>
struct RadiatingBoundary
{
  using Real3 = Vec3<real_t>;

  virtual real_t pulse_s_lower(double t, int d, int p, Real3 x3) = 0;
  virtual real_t pulse_p_lower(double t, int d, int p, Real3 x3) = 0;

  virtual real_t pulse_s_upper(double t, int d, int p, Real3 x3) = 0;
  virtual real_t pulse_p_upper(double t, int d, int p, Real3 x3) = 0;

  // FIXME these are a hack for a specific subclass to work
  virtual void update_cache_lower(double t, int d) {}
  virtual void update_cache_upper(double t, int d) {}
};