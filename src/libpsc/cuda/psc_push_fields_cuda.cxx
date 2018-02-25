
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "push_fields.hxx"

class PushFieldsCuda : PushFieldsBase
{
public:
  // ----------------------------------------------------------------------
  // push_E

  void push_E(struct psc_mfields *mflds_base, double dt_fac) override
  {
    PscMfieldsCuda mf = mflds_base->get_as<PscMfieldsCuda>(JXI, HX + 3);
    
    if (ppsc->domain.gdims[0] == 1) {
      cuda_push_fields_E_yz(mf->cmflds, dt_fac * ppsc->dt);
    } else {
      assert(0);
    }
    
    mf.put_as(mflds_base, EX, EX + 3);
  }
  
  void push_H(struct psc_mfields *mflds_base, double dt_fac) override
  {
    PscMfieldsCuda mf = mflds_base->get_as<PscMfieldsCuda>(JXI, HX + 3);
    
    if (ppsc->domain.gdims[0] == 1) {
      cuda_push_fields_H_yz(mf->cmflds, dt_fac * ppsc->dt);
    } else {
      assert(0);
    }
    
    mf.put_as(mflds_base, HX, HX + 3);
  }
};

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
