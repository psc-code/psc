
#ifndef VPIC_INTERPOLATOR_H
#define VPIC_INTERPOLATOR_H

template<class InterpolatorBase, class FA>
struct VpicInterpolator : InterpolatorBase
{
  typedef InterpolatorBase Base;
  typedef FA FieldArray;

  using Base::Base;
  
  void load(FieldArray& fa)
  {
    TIC ::load_interpolator_array(this, &fa); TOC(load_interpolator, 1);
  }
};

#endif

