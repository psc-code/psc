
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

using fields_single_real_t = float;

struct fields_single_t : fields3d<fields_single_real_t>
{
  using Base = fields3d<fields_single_real_t>;
  using mfields_t = mfields_base<fields_single_t>;

  using Base::Base;
};

struct psc_mfields_single_sub
{
  using fields_t = fields_single_t;
  
  fields_t::real_t **data;
  int ib[3]; //> lower left corner for each patch (incl. ghostpoints)
  int im[3]; //> extent for each patch (incl. ghostpoints)
};

using mfields_single_t = mfields_base<psc_mfields_single_sub>;

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
