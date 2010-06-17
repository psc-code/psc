
#include "psc.h"

#define ALLOC_field_fortran_F77 F77_FUNC_(alloc_field_fortran, ALLOC_FIELD_FORTRAN)
#define SETUP_field_F77 F77_FUNC_(setup_field, SETUP_FIELD)
#define PSC_driver_F77 F77_FUNC_(psc_driver, PSC_DRIVER)

void ALLOC_field_fortran_F77(void);
void SETUP_field_F77(void);
void PSC_driver_F77(void);

void
psc_integrate()
{
  ALLOC_field_fortran_F77();
  SETUP_field_F77();
  PSC_driver_F77();
}
