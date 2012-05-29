
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_fields_as_single.h"

#define psc_bnd_fld_sub_create      psc_bnd_fld_single_create
#define psc_bnd_fld_sub_add_ghosts  psc_bnd_fld_single_add_ghosts
#define psc_bnd_fld_sub_fill_ghosts psc_bnd_fld_single_fill_ghosts

#include "psc_bnd_fld.c"

