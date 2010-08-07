
#include "psc.h"
#include "psc_fields_cuda.h"

void
fields_cuda_get(fields_cuda_t *pf, int mb, int me)
{
  assert(!psc.domain.use_pml);
  int nr_fields = HZ + 1;

  pf->flds = calloc(nr_fields * psc.fld_size, sizeof(*pf->flds));
  
  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_CUDA(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }

  __fields_cuda_get(pf);
}

void
fields_cuda_put(fields_cuda_t *pf, int mb, int me)
{
  __fields_cuda_put(pf);

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_BASE(m, jx,jy,jz) = F3_CUDA(pf, m, jx,jy,jz);
	}
      }
    }
  }

  free(pf->flds);
}

