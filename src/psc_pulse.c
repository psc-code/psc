
#include "psc_pulse_private.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields_private.h"

// ----------------------------------------------------------------------
// psc_pulse_setup

static inline void
_psc_pulse_setup(struct psc_pulse *pulse)
{
  mrc_obj_setup_sub(&pulse->obj);
  pulse->is_setup = true;
}

// ----------------------------------------------------------------------
// forward to subclass

double
psc_pulse_field_s(struct psc_pulse *pulse, double x, double y, double z, double t)
{
  if (!pulse->is_setup) {
    psc_pulse_setup(pulse);
  }
  return psc_pulse_ops(pulse)->field_s(pulse, x, y, z, t);
}

double
psc_pulse_field_p(struct psc_pulse *pulse, double x, double y, double z, double t)
{
  if (!pulse->is_setup) {
    psc_pulse_setup(pulse);
  }
  return psc_pulse_ops(pulse)->field_p(pulse, x, y, z, t);
}

// ======================================================================
// psc_pulse_init

static void
psc_pulse_init()
{
  mrc_class_register_subclass(&mrc_class_psc_pulse, &psc_pulse_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_pulse, &psc_pulse_gauss_ops);
  mrc_class_register_subclass(&mrc_class_psc_pulse, &psc_pulse_flattop_ops);
}

// ======================================================================
// psc_pulse class

#define VAR(x) (void *)offsetof(struct psc_pulse, x)
static struct param psc_pulse_descr[] = {
  { "k"               , VAR(k)            , PARAM_DOUBLE3(0., 0., 0.)       },
  {},
};
#undef VAR

struct mrc_class_psc_pulse mrc_class_psc_pulse = {
  .name             = "psc_pulse",
  .size             = sizeof(struct psc_pulse),
  .param_descr      = psc_pulse_descr,
  .init             = psc_pulse_init,
  .setup            = _psc_pulse_setup,
};

