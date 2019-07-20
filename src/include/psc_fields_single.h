
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

using MfieldsSingle = Mfields<float>;
using MfieldsStateSingle = MfieldsStateFromMfields<MfieldsSingle>;

template<>
struct Mfields_traits<MfieldsSingle>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

template<>
struct Mfields_traits<MfieldsStateSingle>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
