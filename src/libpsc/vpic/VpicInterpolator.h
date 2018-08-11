
#pragma once

// ======================================================================
// VpicInterpolatorOps

template<typename Interpolator, typename FieldArray>
struct VpicInterpolatorOps
{
  static void load(Interpolator& ip, FieldArray& fa)
  {
    TIC ::load_interpolator_array(&ip, &fa); TOC(load_interpolator, 1);
  }
};
