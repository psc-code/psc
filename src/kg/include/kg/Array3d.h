
#ifndef KG_ARRAY3D_H
#define KG_ARRAY3D_H

#include <type_traits>

namespace kg
{

template <bool AOS>
struct Layout
{
  using isAOS = std::integral_constant<bool, AOS>;
};

using LayoutAOS = Layout<true>;
using LayoutSOA = Layout<false>;

} // namespace kg

#endif
