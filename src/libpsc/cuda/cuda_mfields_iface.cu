
#include "cuda_iface.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc_fields_cuda.h"

#if 1
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

psc_mfields_cuda::psc_mfields_cuda(Grid_t& grid, int n_fields, const Int3& ibn)
{
  dprintf("CMFLDS: ctor\n");
  cmflds = new cuda_mfields(grid, n_fields, ibn);
}

psc_mfields_cuda::~psc_mfields_cuda()
{
  dprintf("CMFLDS: dtor\n");
  delete cmflds;
}

fields_single_t psc_mfields_cuda::get_host_fields()
{
  dprintf("CMFLDS: get_host_fields\n");
  return cmflds->get_host_fields();
}

void psc_mfields_cuda::copy_to_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_to_device\n");
  cmflds->copy_to_device(p, h_flds, mb, me);
}

void psc_mfields_cuda::copy_from_device(int p, fields_single_t h_flds, int mb, int me)
{
  dprintf("CMFLDS: copy_from_device\n");
  cmflds->copy_from_device(p, h_flds, mb, me);
}

void psc_mfields_cuda::axpy_comp_yz(int ym, float a, mfields_cuda_t mflds_x, int xm)
{
  dprintf("CMFLDS: axpy_comp_yz\n");
  cmflds->axpy_comp_yz(ym, a, mflds_x->cmflds, xm);
}

void psc_mfields_cuda::zero_comp_yz(int xm)
{
  dprintf("CMFLDS: zero_comp_yz\n");
  cmflds->zero_comp_yz(xm);
}