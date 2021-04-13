
#ifndef FIELDS_HXX
#define FIELDS_HXX

#include "dim.hxx"
#include "kg/Vec3.h"

#include <gtensor/gtensor.h>

// ======================================================================
// Fields3d
//
// wrappers a gtensor expression, shifting the offset from zero and
// setting to 0 indices in the invariant direction

template <typename F, typename D = dim_xyz>
class Fields3d
{
public:
  using fields_t = F;
  using value_type = typename fields_t::value_type;
  using shape_type = typename fields_t::shape_type;
  using dim = D;

  GT_INLINE Fields3d(const fields_t& e, const Int3& ib) : e_(e), ib_(ib) {}

  GT_INLINE shape_type shape() const { return e_.shape(); }
  GT_INLINE int shape(int d) const { return e_.shape(d); }
  GT_INLINE Int3 ib() const { return ib_; }

  GT_INLINE const value_type& operator()(int m, int _i, int _j, int _k) const
  {
    int i = dim::InvarX::value ? 0 : _i - ib_[0];
    int j = dim::InvarY::value ? 0 : _j - ib_[1];
    int k = dim::InvarZ::value ? 0 : _k - ib_[2];

    return e_(i, j, k, m);
  }

  GT_INLINE value_type& operator()(int m, int _i, int _j, int _k)
  {
    int i = dim::InvarX::value ? 0 : _i - ib_[0];
    int j = dim::InvarY::value ? 0 : _j - ib_[1];
    int k = dim::InvarZ::value ? 0 : _k - ib_[2];

    return e_(i, j, k, m);
  }

private:
  fields_t e_;
  Int3 ib_;
};

template <typename D, typename F>
GT_INLINE auto make_Fields3d(const F& f)
{
  return Fields3d<typename F::Storage, D>(f.storage(), f.ib());
}

#endif
