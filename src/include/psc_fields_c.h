
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

struct fields_c_t : fields3d<double>
{
  using Base = fields3d<double>;
  using mfields_t = mfields3d<fields_c_t>;

  using Base::Base;
};

using mfields_c_t = mfields3d<fields_c_t>;

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
