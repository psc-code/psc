
#ifndef KG_ARRAY3D_H
#define KG_ARRAY3D_H

#include <kg/Vec3.h>

#include <type_traits>

namespace kg
{

struct LayoutAOS
{
  using isAOS = std::true_type;
};

struct LayoutSOA
{
  using isAOS = std::false_type;
};

// ======================================================================
// Box3

class Box3
{
public:
  Box3(const Int3& ib, const Int3& im)
    : ib_{ib}, im_{im}
  {}

  KG_INLINE int size() const { return im_[0] * im_[1] * im_[2]; }

  KG_INLINE const Int3& ib() const { return ib_; }
  KG_INLINE const Int3& im() const { return im_; }
  KG_INLINE int ib(int d) const { return ib_[d]; }
  KG_INLINE int im(int d) const { return im_[d]; }

private:
  Int3 ib_;
  Int3 im_;
};

} // namespace kg

#endif
