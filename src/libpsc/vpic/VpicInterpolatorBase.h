
#ifndef VPIC_INTERPOLATOR_BASE_H
#define VPIC_INTERPOLATOR_BASE_H

// ======================================================================
// VpicInterpolatorBase

template<class G>
struct VpicInterpolatorBase : interpolator_array_t {
  typedef G Grid;
  typedef interpolator_t Element;

  static VpicInterpolatorBase* create(Grid *grid)
  {
    return static_cast<VpicInterpolatorBase*>(new_interpolator_array(grid));
  }
  
  static void destroy(VpicInterpolatorBase* interpolator)
  {
    delete_interpolator_array(interpolator);
  }
  
  Element  operator[](int idx) const { return i[idx]; }
  Element& operator[](int idx)       { return i[idx]; }

  Element* data() { return i; }
};

#endif

