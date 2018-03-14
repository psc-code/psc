
#include "cuda_iface.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc_fields_cuda.h"

#if 1
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

MfieldsCuda::MfieldsCuda(const Grid_t& grid, int n_fields, const Int3& ibn)
  : MfieldsBase(grid, n_fields)
{
  dprintf("CMFLDS: ctor\n");
  cmflds = new cuda_mfields(grid, n_fields, ibn);
}

MfieldsCuda::~MfieldsCuda()
{
  dprintf("CMFLDS: dtor\n");
  delete cmflds;
}

fields_single_t MfieldsCuda::get_host_fields()
{
  dprintf("CMFLDS: get_host_fields\n");
  return cmflds->get_host_fields();
}

void MfieldsCuda::copy_to_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_to_device\n");
  cmflds->copy_to_device(p, h_flds, mb, me);
}

void MfieldsCuda::copy_from_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_from_device\n");
  cmflds->copy_from_device(p, h_flds, mb, me);
}

void MfieldsCuda::axpy_comp_yz(int ym, float a, PscMfieldsCuda mflds_x, int xm)
{
  dprintf("CMFLDS: axpy_comp_yz\n");
  cmflds->axpy_comp_yz(ym, a, mflds_x->cmflds, xm);
}

void MfieldsCuda::zero_comp(int m)
{
  dprintf("CMFLDS: zero_comp\n");
  assert(grid_->gdims[0] == 1);
  cmflds->zero_comp_yz(m);
}

void MfieldsCuda::zero()
{
  dprintf("CMFLDS: zero\n");
  for (int m = 0; m < cmflds->n_fields; m++) {
    zero_comp(m);
  }
}