
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#define FTYPE FTYPE_SINGLE
#include "psc_fields_common.h"
#undef FTYPE

#include "fields_traits.hxx"

template<>
struct fields_traits<fields_single_t>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
