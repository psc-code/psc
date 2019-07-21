
#ifndef KG_ARRAY3D_H
#define KG_ARRAY3D_H

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

} // namespace kg

#endif
