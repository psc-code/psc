
#include "push_fields.hxx"
#include "dim.hxx"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

struct PushFieldsCuda : PushFieldsBase
{
  void push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag);
  void push_E(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag);

  void push_H(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag);
  void push_H(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag);
};
