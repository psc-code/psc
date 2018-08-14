
#ifndef VPIC_INTERPOLATOR_BASE_H
#define VPIC_INTERPOLATOR_BASE_H

#include "sf_interface/sf_interface.h"

// ======================================================================
// VpicInterpolatorBase

template<class G>
struct VpicInterpolatorBase : interpolator_array_t {
  typedef G Grid;
  typedef interpolator_t Element;

  Element  operator[](int idx) const { return i[idx]; }
  Element& operator[](int idx)       { return i[idx]; }

  Element* data() { return i; }

  Grid* grid() { return static_cast<Grid*>(g); }
};

#endif

