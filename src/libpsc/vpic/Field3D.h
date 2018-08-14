
#ifndef FIELD3D_H
#define FIELD3D_H

#include "psc_vpic_bits.h"

// ======================================================================
// Field3D
//
// A class to accelerate 3-d field access
// It works for VpicFieldArray and VpicInterpolatorArray, but should
// be generalized to work without knowing their internals

template<class FA>
struct Field3D {
  typedef FA Array;
  typedef typename FA::Element Element;
  
  Field3D(Array& fa)
    : sx_(fa.grid()->nx + 2), sy_(fa.grid()->ny + 2),
      f_(fa.data())
  {
  }

  int voxel(int i, int j, int k) const
  {
    return i + sx_ * (j + sy_ * (k));
  }

  // FIXME, this is an odd mix of two interfaces
  // first field3d just being used for grid information,
  // and providing access to the whole struct
  // This interface can't easily be converted to SOA
  Element& operator()(Array &fa, int i, int j, int k)
  {
    return fa.data()[voxel(i,j,k)];
  }
  
  Element operator()(Array &fa, int i, int j, int k) const
  {
    return fa.data()[voxel(i,j,k)];
  }

  Element& operator()(int i, int j, int k)
  {
    return f_[voxel(i,j,k)];
  }

  Element operator()(int i, int j, int k) const
  {
    return f_[voxel(i,j,k)];
  }

  // second, access to a particular component, but this one is
  // for the specific Array used at construction time
  float& operator()(int m, int i, int j, int k)
  {
    float * RESTRICT f = reinterpret_cast<float *>(f_);
    return f[m + Array::N_COMP * voxel(i,j,k)];
  }
  
  float operator()(int m, int i, int j, int k) const
  {
    float * RESTRICT f = reinterpret_cast<float *>(f_);
    return f[m + Array::N_COMP * voxel(i,j,k)];
  }

private:
  int sx_, sy_;
  Element * RESTRICT f_;
};

#endif
