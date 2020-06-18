
#pragma once

#include <kg/io.h>

#include "FileAdios2.h"

namespace kg
{
namespace io
{

// ======================================================================
// IOAdios2

class IOAdios2
{
public:
  IOAdios2();
  IOAdios2(const std::string& config);

  File openFile(const std::string& name, const Mode mode,
                MPI_Comm comm = MPI_COMM_WORLD,
                const std::string& io_name = {});
  Engine open(const std::string& name, const Mode mode,
              MPI_Comm comm = MPI_COMM_WORLD, const std::string& io_name = {});

private:
  adios2::ADIOS ad_;
};

} // namespace io
} // namespace kg

#include "IOAdios2.inl"
