
#include "psc_push_fields_impl.hxx"

#include "psc_push_fields_iface.h"

#include "psc.h"

void psc_push_fields_single_push_E_xz(struct psc_push_fields* push, fields_single_t flds,
				      struct psc* psc, double dt_fac)
{
  using Fields = Fields3d<fields_single_real_t, fields_single_t, DIM_XZ>;
  psc_push_fields_push_E<Fields>(push, flds, psc, dt_fac);
}

void psc_push_fields_single_push_H_xz(struct psc_push_fields* push, fields_single_t flds,
				      struct psc* psc, double dt_fac)
{
  using Fields = Fields3d<fields_single_real_t, fields_single_t, DIM_XZ>;
  psc_push_fields_push_H<Fields>(push, flds, psc, dt_fac);
}

