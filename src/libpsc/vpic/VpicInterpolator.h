
#ifndef VPIC_INTERPOLATOR_H
#define VPIC_INTERPOLATOR_H

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

template<class InterpolatorBase, class FieldArray>
struct VpicInterpolator : InterpolatorBase
{
  using Base = InterpolatorBase;
  using Self = VpicInterpolator<InterpolatorBase, FieldArray>;

  using Base::Base;

  void load(FieldArray& fa) { VpicInterpolatorOps<Self, FieldArray>::load(*this, fa); }
};

#endif

