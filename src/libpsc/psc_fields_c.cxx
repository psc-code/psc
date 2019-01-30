
#include "psc.h"
#include "psc_fields_c.h"

#define PFX(x) psc_fields_c_ ## x
#define MPFX(x) psc_mfields_c_ ## x
#define MFIELDS MfieldsC

template<> const MfieldsBase::Convert MfieldsC::convert_to_{};
template<> const MfieldsBase::Convert MfieldsC::convert_from_{};

template<> const MfieldsStateBase::Convert MfieldsStateDouble::convert_to_{};
template<> const MfieldsStateBase::Convert MfieldsStateDouble::convert_from_{};

#include "psc_fields_common.cxx"

