
#include "psc_pulse_private.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields_private.h"

// ----------------------------------------------------------------------
// psc_p_pulse_x1

real
psc_p_pulse_x1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_x1)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_x1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_x1

real
psc_s_pulse_x1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_x1)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_x1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_x2

real
psc_p_pulse_x2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_x2)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_x2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_x2

real
psc_s_pulse_x2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_x2)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_x2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y1

real
psc_p_pulse_y1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_y1)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_y1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_y1

real
psc_s_pulse_y1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_y1)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_y1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y2

real
psc_p_pulse_y2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_y2)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_y2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_y2

real
psc_s_pulse_y2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_y2)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_y2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z1

real
psc_p_pulse_z1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_z1)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_z1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z1

real
psc_s_pulse_z1(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_z1)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_z1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z2

real
psc_p_pulse_z2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_z2)
    return 0.;

  return psc_pulse_field_p(bnd_fields->pulse_z2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z2

real
psc_s_pulse_z2(real x, real y, real z, real t)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  if (!bnd_fields->pulse_z2)
    return 0.;

  return psc_pulse_field_s(bnd_fields->pulse_z2, x, y, z, t);
}


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

