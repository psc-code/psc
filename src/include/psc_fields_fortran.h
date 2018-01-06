
#ifndef PSC_FIELD_FORTRAN_H
#define PSC_FIELD_FORTRAN_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_FORTRAN
#include "psc_fields_common.h"
#undef FTYPE

template<>
inline fields_fortran_t mfields_fortran_t::operator[](int p)
{
  fields_fortran_t psc_mfields_fortran_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_fortran_get_field_t(mflds_, p);
}

#endif
