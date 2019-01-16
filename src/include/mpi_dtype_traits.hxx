
#pragma once

#include <mpi.h>

template<typename T>
struct MpiDtypeTraits;

template<>
struct MpiDtypeTraits<float>
{
  static MPI_Datatype value() { return MPI_FLOAT; }
};

template<>
struct MpiDtypeTraits<double>
{
  static MPI_Datatype value() { return MPI_DOUBLE; }
};
