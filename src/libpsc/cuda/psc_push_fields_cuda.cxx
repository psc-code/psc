
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "push_fields.hxx"
#include "dim.hxx"

struct PushFieldsCuda : PushFieldsBase
{
  // ----------------------------------------------------------------------
  // push_E

  void push_E(MfieldsCuda& mflds, double dt_fac, dim_yz tag)
  {
    cuda_push_fields_E_yz(mflds.cmflds, dt_fac * ppsc->dt);
  }

  // dispatch -- FIXME, mostly same as non-cuda dispatch
  void push_E(PscMfieldsBase mflds_base, double dt_fac) override
  {
    auto& mflds = mflds_base->get_as<MfieldsCuda>(JXI, HX + 3);
    
    const auto& grid = mflds.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    if (invar == Bool3{true, false, false}) {
      push_E(mflds, dt_fac, dim_yz{});
    } else {
      assert(0);
    }
    
    mflds_base->put_as(mflds, EX, EX + 3);
  }

  // ----------------------------------------------------------------------
  // push_H

  void push_H(MfieldsCuda& mflds, double dt_fac, dim_yz tag)
  {
    cuda_push_fields_H_yz(mflds.cmflds, dt_fac * ppsc->dt);
  }

  // dispatch -- FIXME, mostly same as non-cuda dispatch
  void push_H(PscMfieldsBase mflds_base, double dt_fac) override
  {
    auto& mflds = mflds_base->get_as<MfieldsCuda>(JXI, HX + 3);
    
    const auto& grid = mflds.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    if (invar == Bool3{true, false, false}) {
      push_H(mflds, dt_fac, dim_yz{});
    } else {
      assert(0);
    }

    mflds_base->put_as(mflds, HX, HX + 3);
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
