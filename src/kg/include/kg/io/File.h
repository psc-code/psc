
#pragma once

#include "FileBase.h"

namespace kg
{
namespace io
{

// ======================================================================
// File

class File
{
public:
  explicit File(FileBase* impl);
  File(File&& other) = default;
  ~File();
  
  void close();
  void performPuts();
  void performGets();

  template <typename T>
  void putVariable(const std::string& name, const T* data, Mode launch,
                   const Dims& shape, const Extents& selection,
                   const Extents& memory_selection);

  template <typename T>
  void getVariable(const std::string& name, T* data, Mode launch,
                   const Extents& selection,
                   const Extents& memory_selection);

  Dims shapeVariable(const std::string& name) const;

  template <typename T>
  void getAttribute(const std::string& name, T* data);

  template <typename T>
  void putAttribute(const std::string& name, const T* data, size_t size);

  size_t sizeAttribute(const std::string& name) const;

private:
  std::unique_ptr<FileBase> impl_;
};

} // namespace io
} // namespace kg

#include "File.inl"

