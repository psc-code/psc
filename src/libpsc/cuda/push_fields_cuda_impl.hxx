
#include "push_fields.hxx"
#include "dim.hxx"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

struct PushFieldsCuda : PushFieldsBase
{
  // ----------------------------------------------------------------------
  // push_E

  void push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
  {
    cuda_push_fields_E_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
  }

  void push_E(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
  {
    cuda_push_fields_E_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
  }

  // ----------------------------------------------------------------------
  // push_H

  void push_H(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
  {
    cuda_push_fields_H_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
  }

  void push_H(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
  {
    cuda_push_fields_H_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
  }
};
