
#pragma once

#include "mpark/variant.hpp"

namespace kg
{
namespace io
{

// ======================================================================
// FileBase

class FileBase
{
public:
  using TypePointer =
    mpark::variant<int*, unsigned long*, double*, std::string*>;
  using TypeConstPointer = mpark::variant<const int*, const unsigned long*,
                                          const double*, const std::string*>;

  virtual ~FileBase() = default;

  virtual void performPuts() = 0;
  virtual void performGets() = 0;

  virtual void putVariable(const std::string& name, TypeConstPointer data,
                           Mode launch, const Dims& shape,
                           const Extents& selection,
                           const Extents& memory_selection) = 0;

  virtual void getVariable(const std::string& name, TypePointer data,
                           Mode launch, const Extents& selection,
                           const Extents& memory_selection) = 0;
  virtual Dims shapeVariable(const std::string& name) const = 0;

  virtual void getAttribute(const std::string& name, TypePointer data) = 0;
  virtual void putAttribute(const std::string& name, TypeConstPointer data,
                            size_t size) = 0;
  virtual size_t sizeAttribute(const std::string& name) const = 0;
};

}
} // namespace kg
