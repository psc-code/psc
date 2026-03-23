#include "kg/Vec3.h"

template <typename real_t>
class RadiatingBoundary
{
  virtual real_t pulse_s_lower(int d, int p, Int3 i3) = 0;
  virtual real_t pulse_p_lower(int d, int p, Int3 i3) = 0;

  virtual real_t pulse_s_upper(int d, int p, Int3 i3) = 0;
  virtual real_t pulse_p_upper(int d, int p, Int3 i3) = 0;
};