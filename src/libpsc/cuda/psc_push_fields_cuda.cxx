
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "push_fields_cuda_impl.hxx"

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops_cuda : psc_push_fields_ops {
  using Wrapper = PscPushFieldsWrapper<PushFieldsCuda>;
  psc_push_fields_ops_cuda() {
    name                  = "cuda";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_push_fields_cuda_ops;

