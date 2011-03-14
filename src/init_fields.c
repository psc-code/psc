
#include "psc.h"
#include "psc_case_private.h"

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
  if (_psc_case->Case) {
    psc_case_init_field(_psc_case->Case, flds);
    if (psc.domain.use_pml) {
      psc_init_field_pml(flds);
    }
  } else {
    INIT_field();
  }
}

// ======================================================================
// Fortran glue

#define C_INIT_field_F77 F77_FUNC_(c_init_field, C_INIT_FIELD)
#define INIT_field_F77 F77_FUNC_(init_field, INIT_FIELD)

void INIT_field_F77(void);

void
INIT_field()
{
  INIT_field_F77();
}

void
C_INIT_field_F77()
{
  psc_init_field(&psc.flds);
}
