
#pragma once

#include <kg/io.h>

namespace
{

inline kg::io::Dims makeDims(int m, const Int3& dims)
{
  return kg::io::Dims{static_cast<size_t>(m), static_cast<size_t>(dims[2]),
                      static_cast<size_t>(dims[1]),
                      static_cast<size_t>(dims[0])};
}

} // namespace
