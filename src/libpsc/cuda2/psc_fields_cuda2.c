
#include "psc_fields_cuda2.h"

// ======================================================================
// psc_fields: subclass "cuda2"
  
struct psc_fields_ops psc_fields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_fields_cuda2),
#if 0
  .methods               = psc_fields_cuda_methods,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_cuda_write,
#endif
  .axpy_comp             = psc_fields_cuda_axpy_comp,
  .zero_comp             = psc_fields_cuda_zero_comp,
#endif
};

// ======================================================================
// psc_mfields: subclass "cuda2"
  
struct psc_mfields_ops psc_mfields_cuda2_ops = {
  .name                  = "cuda2",
  .size                  = sizeof(struct psc_mfields_cuda2),
#if 0
  .setup                 = psc_mfields_cuda2_setup,
  .destroy               = psc_mfields_cuda2_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_mfields_cuda2_read,
#endif
#endif
};
