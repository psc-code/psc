
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_C
#include "psc_fields_common.h"
#undef FTYPE

#include "fields_traits.hxx"

template<>
inline fields_c_t mfields_c_t::operator[](int p)
{
  fields_c_t psc_mfields_c_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_c_get_field_t(mflds_, p);
}

template<>
struct fields_traits<fields_c_t>
{
  static constexpr const char* name = "c";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
