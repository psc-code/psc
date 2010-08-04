
#include "psc.h"
#include "psc_fields_fortran.h"

void
psc_fields_fortran_alloc(psc_fields_fortran_t *pf)
{
  f_real **fields = ALLOC_field();
  for (int i = 0; i < NR_FIELDS; i++) {
    pf->flds[i] = fields[i];
  }
}

void
psc_fields_fortran_free(psc_fields_fortran_t *pf)
{
  FREE_field();
  for (int i = 0; i < NR_FIELDS; i++) {
    pf->flds[i] = NULL;
  }
}

void
psc_fields_fortran_zero(psc_fields_fortran_t *pf, int m)
{
  memset(pf->flds[m], 0, psc.fld_size * sizeof(f_real));
}

