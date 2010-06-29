
#include "psc.h"

#include <string.h>

void INIT_field(void);

void
init_field()
{
  if (psc.case_ops) {
    if (psc.case_ops->init_field) {
      psc.case_ops->init_field();
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
  init_field();
}
