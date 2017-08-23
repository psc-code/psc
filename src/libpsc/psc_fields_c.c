
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

static struct mrc_obj_method psc_fields_c_methods[] = {
  {}
};

#include "psc_fields_common.c"

