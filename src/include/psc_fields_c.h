
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_C
#include "psc_fields_common.h"
#undef FTYPE

#include "fields_traits.hxx"

template<>
struct fields_traits<fields_c_t>
{
  static constexpr const char* name = "c";
};

#endif
