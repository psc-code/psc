
#pragma once

// ======================================================================
// VpicInterpolatorOps

template<typename Interpolator, typename FieldArray>
struct VpicInterpolatorOps
{
  static void load(Interpolator& ip, FieldArray& fa) { ::load_interpolator_array(&ip, &fa); }
};
