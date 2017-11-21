
#ifndef VPIC_INTERPOLATOR_OPS_H
#define VPIC_INTERPOLATOR_OPS_H

template<class I, class F>
struct VpicInterpolatorOps {
  typedef I Interpolator;
  typedef F FieldArray;

  void load_interpolator_array(Interpolator *interpolator,
			       FieldArray *vmflds)
  {
    TIC ::load_interpolator_array(interpolator, vmflds); TOC(load_interpolator, 1);
  }
};

#endif

