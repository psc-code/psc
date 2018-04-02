
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "push_fields_cuda_impl.hxx"

static void
psc_push_fields_cuda_setup(struct psc_push_fields *push)
{
  PscPushFields<PushFieldsCuda> pushf(push);
  new(pushf.sub()) PushFieldsCuda;
}

static void
psc_push_fields_cuda_destroy(struct psc_push_fields *push)
{
  PscPushFields<PushFieldsCuda> pushf(push);
  pushf.sub()->~PushFieldsCuda();
}

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops_cuda : psc_push_fields_ops {
  psc_push_fields_ops_cuda() {
    name                  = "cuda";
    size                  = sizeof(PushFieldsCuda),
    setup                 = psc_push_fields_cuda_setup;
    destroy               = psc_push_fields_cuda_destroy;
  }
} psc_push_fields_cuda_ops;
