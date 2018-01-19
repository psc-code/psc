
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

using fields_c_real_t = double;

struct fields_c_t : fields3d<fields_c_real_t>
{
  using Base = fields3d<fields_c_real_t>;
  using mfields_t = mfields_base<fields_c_t>;

  using Base::Base;
};

struct psc_mfields_c_sub
{
  using fields_t = fields_c_t;
  
  fields_t::real_t **data;
  int ib[3]; //> lower left corner for each patch (incl. ghostpoints)
  int im[3]; //> extent for each patch (incl. ghostpoints)
};

using mfields_c_t = mfields_base<psc_mfields_c_sub>;

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
