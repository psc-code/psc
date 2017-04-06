
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_fld_as_double.h>

#define MT MT_FORMULATION_SCONS
#define SFX(x) x ## _scons

#include "ggcm_mhd_primitive_common.c"
