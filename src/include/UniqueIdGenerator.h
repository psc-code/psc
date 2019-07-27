
#pragma once

#include "particle.h"

#include <mpi.h>

namespace psc
{
namespace particle
{

class UniqueIdGenerator
{
public:
  UniqueIdGenerator(MPI_Comm comm)
  {
    MPI_Comm_rank(comm, &rank_);
    MPI_Comm_size(comm, &size_);
  }

  Id operator()() { return id_++ * size_ + rank_; }

private:
  int rank_;
  int size_;
  Id id_ = 0;
};

} // namespace particle
} // namespace psc
