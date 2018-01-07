
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

struct fields_single_t : fields3d<float>
{
  using Base = fields3d<float>;
  using mfields_t = mfields3d<fields_single_t>;

  using Base::Base;
};

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
