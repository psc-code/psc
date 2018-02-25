
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

using fields_c_real_t = double;

struct fields_c_t : fields3d<fields_c_real_t>
{
  using Base = fields3d<fields_c_real_t>;
  using mfields_t = PscMfields<fields_c_t>;

  using Base::Base;
};

using psc_mfields_c_sub = psc_mfields_<fields_c_t>;
using PscMfieldsC = PscMfields<psc_mfields_c_sub>;

template<>
struct fields_traits<fields_c_t>
{
  static constexpr const char* name = "c";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
