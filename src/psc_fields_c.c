
#include "psc.h"
#include "psc_fields_c.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// FIXME, we're allocating all fields even if the PML ones aren't
// going to be used

void
fields_c_alloc(fields_c_t *pf, int ib[3], int ie[3], int nr_comp)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->flds = calloc(nr_comp * psc.fld_size, sizeof(*pf->flds));
  pf->name = calloc(nr_comp, sizeof(*pf->name));
  for (int m = 0; m < nr_comp; m++) {
    pf->name[m] = NULL;
  }
}

void
fields_c_free(fields_c_t *pf)
{
  free(pf->flds);
  for (int m = 0; m < pf->nr_comp; m++) {
    free(pf->name[m]);
  }
  free(pf->name);
}

#if FIELDS_BASE == FIELDS_C

void
fields_c_get(fields_c_t *pf, int mb, int me)
{
  *pf = psc.pf;
}

void
fields_c_put(fields_c_t *pf, int mb, int me)
{
  pf->flds = NULL;
}

#else

void
fields_c_get(fields_c_t *pf, int mb, int me)
{
  fields_c_alloc(pf, psc.ilg, psc.ihg, NR_FIELDS);

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
fields_c_put(fields_c_t *pf, int mb, int me)
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

  fields_c_free(pf);
  pf->flds = NULL;
}

#endif

void
fields_c_zero(fields_c_t *pf, int m)
{
  memset(&F3_C(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_c_real_t));
}

void
fields_c_set(fields_c_t *pf, int m, fields_c_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_C(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_c_copy(fields_c_t *pf, int m_to, int m_from)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_C(pf, m_to, jx, jy, jz) = F3_C(pf, m_from, jx, jy, jz);
      }
    }
  }
}
