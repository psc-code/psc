
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

struct fields_single_t : fields3d<float>
{
  using Base = fields3d<float>;

  using Base::Base;
};

using MfieldsSingle = Mfields<fields_single_t>;
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
