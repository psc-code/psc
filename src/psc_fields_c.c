
#include "psc.h"
#include "psc_fields_c.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

void
psc_fields_c_get(psc_fields_c_t *pf, int mb, int me)
{
  pf->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));

  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_C(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }
}

void
psc_fields_c_put(psc_fields_c_t *pf, int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_BASE(m, jx,jy,jz) = F3_C(pf, m, jx,jy,jz);
	}
      }
    }
  }

  free(pf->flds);
}

