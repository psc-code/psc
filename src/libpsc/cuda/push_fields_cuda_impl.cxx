
#include "push_fields_cuda_impl.hxx"

// ----------------------------------------------------------------------
// push_E

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  cuda_push_fields_E_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_E_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

// ----------------------------------------------------------------------
// push_H

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  cuda_push_fields_H_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_H_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}
