
#include "psc.h"
#include "psc_fields_cuda.h"

void
fields_cuda_get(fields_cuda_t *pf, int mb, int me)
{
  assert(!psc.domain.use_pml);
  int nr_fields = HZ + 1;

  pf->flds = calloc(nr_fields * psc.fld_size, sizeof(*pf->flds));
  
  for (int m = mb; m < me; m++) {
    foreach_3d_g(jx, jy, jz) {
      F3_CUDA(pf, m, jx,jy,jz) = XF3_BASE(&psc.pf, m, jx,jy,jz);
    } foreach_3d_g_end;
  }

  __fields_cuda_get(pf);
}

void
fields_cuda_put(fields_cuda_t *pf, int mb, int me)
{
  __fields_cuda_put(pf);

  for (int m = mb; m < me; m++) {
    foreach_3d_g(jx, jy, jz) {
      XF3_BASE(&psc.pf, m, jx,jy,jz) = F3_CUDA(pf, m, jx,jy,jz);
    } foreach_3d_g_end;
  }

  free(pf->flds);
}

