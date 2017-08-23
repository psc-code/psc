
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "psc_fields_inc.h"

#include "psc_fields_common.c"

// ======================================================================
// psc_fields: subclass "c"
  
struct psc_fields_ops psc_fields_c_ops = {
  .name                  = "c",
  .setup                 = psc_fields_c_setup,
  .destroy               = psc_fields_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_fields_c_read,
  .write                 = psc_fields_c_write,
#endif
  .zero_comp             = psc_fields_c_zero_comp,
  .set_comp              = psc_fields_c_set_comp,
  .scale_comp            = psc_fields_c_scale_comp,
  .copy_comp             = psc_fields_c_copy_comp,
  .axpy_comp             = psc_fields_c_axpy_comp,
};

