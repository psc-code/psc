
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define PFX(x) psc_fields_c_ ## x
#define MPFX(x) psc_mfields_c_ ## x
#define MFIELDS MfieldsC

template<> const MfieldsBase::Convert MfieldsC::convert_to_ = {};
template<> const MfieldsBase::Convert MfieldsC::convert_from_ = {};

#include "psc_fields_common.cxx"

