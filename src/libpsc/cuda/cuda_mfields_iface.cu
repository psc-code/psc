
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

