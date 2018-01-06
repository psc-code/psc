
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include "fields_traits.hxx"

#define FTYPE FTYPE_SINGLE
#include "psc_fields_common.h"
#undef FTYPE

using mfields_single_t = mfields3d<fields_single_t>;

template<>
inline fields_single_t mfields_single_t::operator[](int p)
{
  fields_single_t psc_mfields_single_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_single_get_field_t(mflds_, p);
}

template<>
struct fields_traits<fields_single_t>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
