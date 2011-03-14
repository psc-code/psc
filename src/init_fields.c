
#include "psc.h"

#include <string.h>

void INIT_field(void);

static void
psc_init_field_pml(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_copy(&flds->f[p], DX, EX);
    fields_base_copy(&flds->f[p], DY, EY);
    fields_base_copy(&flds->f[p], DZ, EZ);
    fields_base_copy(&flds->f[p], BX, HX);
    fields_base_copy(&flds->f[p], BY, HY);
    fields_base_copy(&flds->f[p], BZ, HZ);
    fields_base_set(&flds->f[p], EPS, 1.);
    fields_base_set(&flds->f[p], MU, 1.);
  }
}

void
psc_init_field(mfields_base_t *flds)
{
  psc_case_init_field(_psc_case, flds);
  if (psc.domain.use_pml) {
    psc_init_field_pml(flds);
  }
}

