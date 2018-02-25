
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

struct fields_c_t : fields3d<double>
{
  using Base = fields3d<double>;

  using Base::Base;
};

using MfieldsC = Mfields<fields_c_t>;
using PscMfieldsC = PscMfields<MfieldsC>;

template<>
struct fields_traits<fields_c_t>
{
  static constexpr const char* name = "c";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
