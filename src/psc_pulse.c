
#include "psc.h"
#include "psc_glue.h"

// ----------------------------------------------------------------------
// psc_p_pulse_x1

real
psc_p_pulse_x1(real x, real y, real z, real t)
{
  if (!psc.pulse_x1)
    return 0.;

  return psc_pulse_field_p(psc.pulse_x1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_x1

real
psc_s_pulse_x1(real x, real y, real z, real t)
{
  if (!psc.pulse_x1)
    return 0.;

  return psc_pulse_field_s(psc.pulse_x1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_x2

real
psc_p_pulse_x2(real x, real y, real z, real t)
{
  if (!psc.pulse_x2)
    return 0.;

  return psc_pulse_field_p(psc.pulse_x2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_x2

real
psc_s_pulse_x2(real x, real y, real z, real t)
{
  if (!psc.pulse_x2)
    return 0.;

  return psc_pulse_field_s(psc.pulse_x2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y1

real
psc_p_pulse_y1(real x, real y, real z, real t)
{
  if (!psc.pulse_y1)
    return 0.;

  return psc_pulse_field_p(psc.pulse_y1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_y1

real
psc_s_pulse_y1(real x, real y, real z, real t)
{
  if (!psc.pulse_y1)
    return 0.;

  return psc_pulse_field_s(psc.pulse_y1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y2

real
psc_p_pulse_y2(real x, real y, real z, real t)
{
  if (!psc.pulse_y2)
    return 0.;

  return psc_pulse_field_p(psc.pulse_y2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_y2

real
psc_s_pulse_y2(real x, real y, real z, real t)
{
  if (!psc.pulse_y2)
    return 0.;

  return psc_pulse_field_s(psc.pulse_y2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z1

real
psc_p_pulse_z1(real x, real y, real z, real t)
{
  if (!psc.pulse_z1)
    return 0.;

  return psc_pulse_field_p(psc.pulse_z1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z1

real
psc_s_pulse_z1(real x, real y, real z, real t)
{
  if (!psc.pulse_z1)
    return 0.;

  return psc_pulse_field_s(psc.pulse_z1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z2

real
psc_p_pulse_z2(real x, real y, real z, real t)
{
  if (!psc.pulse_z2)
    return 0.;

  return psc_pulse_field_p(psc.pulse_z2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z2

real
psc_s_pulse_z2(real x, real y, real z, real t)
{
  if (!psc.pulse_z2)
    return 0.;

  return psc_pulse_field_s(psc.pulse_z2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_pulse_setup

static inline void
_psc_pulse_setup(struct psc_pulse *pulse)
{
  if (pulse->ops->setup) {
    pulse->ops->setup(pulse);
  }
  pulse->is_setup = true;
}

// ----------------------------------------------------------------------
// psc_pulse_ini

void
psc_pulse_ini(struct psc_pulse *pulse, struct psc_pulse_ops *ops, void *prm)
{
  pulse->ops = ops;
  
  if (ops->ctx_size) {
    pulse->ctx = malloc(ops->ctx_size);
    memset(pulse->ctx, 0, ops->ctx_size);
  }

  struct param *descr = ops->ctx_descr;
  if (descr) {
    char name[strlen(ops->name) + 7];
    sprintf(name, "Pulse %s", ops->name);
    if (prm) { // custom defaults were passed
      memcpy(pulse->ctx, prm, ops->ctx_size);
      mrc_params_parse_nodefault(pulse->ctx, descr, name, MPI_COMM_WORLD);
    } else {
      mrc_params_parse(pulse->ctx, descr, name, MPI_COMM_WORLD);
    }
    mrc_params_print(pulse->ctx, descr, name, MPI_COMM_WORLD);
  }
}

// ======================================================================
// psc_pulse_init

static void
psc_pulse_init()
{
  //  mrc_class_register_subclass(&mrc_class_psc_pulse, &psc_pulse_gauss_ops);
}

// ======================================================================
// psc_pulse class

struct mrc_class_psc_pulse mrc_class_psc_pulse = {
  .name             = "psc_pulse",
  .size             = sizeof(struct psc_pulse),
  .init             = psc_pulse_init,
  .setup            = _psc_pulse_setup,
};

