
#ifndef PSC_FIELD_FORTRAN_H
#define PSC_FIELD_FORTRAN_H

#include "fields3d.hxx"

struct fields_fortran_t : fields3d<double>
{
  using Base = fields3d<double>;
  using mfields_t = mfields3d<fields_fortran_t>;

  using Base::Base;
};

using mfields_fortran_t = mfields3d<fields_fortran_t>;

template<>
inline fields_fortran_t mfields_fortran_t::operator[](int p)
{
  fields_fortran_t psc_mfields_fortran_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_fortran_get_field_t(mflds_, p);
}

#endif
