
#include "psc.h"
#include "psc_glue.h"

// ----------------------------------------------------------------------
// psc_p_pulse_x1

real
psc_p_pulse_x1(real x, real y, real z, real t)
{
  if (!psc.pulse_p_x1)
    return 0.;

  return psc_pulse_field(psc.pulse_p_x1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_x2

real
psc_p_pulse_x2(real x, real y, real z, real t)
{
  if (!psc.pulse_p_x2)
    return 0.;

  return psc_pulse_field(psc.pulse_p_x2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y1

real
psc_p_pulse_y1(real x, real y, real z, real t)
{
  if (!psc.pulse_p_y1)
    return 0.;

  return psc_pulse_field(psc.pulse_p_y1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_y2

real
psc_p_pulse_y2(real x, real y, real z, real t)
{
  if (!psc.pulse_p_y2)
    return 0.;

  return psc_pulse_field(psc.pulse_p_y2, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z1

real
psc_p_pulse_z1(real x, real y, real z, real t)
{
  if (!psc.pulse_p_z1)
    return 0.;

  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_p_z1) { // default to Fortran
    return PSC_p_pulse_z1(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_p_z1, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z1

real
psc_s_pulse_z1(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_s_z1) { // default to Fortran
    return PSC_s_pulse_z1(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_s_z1, x, y, z, t);
}

// ----------------------------------------------------------------------
// psc_p_pulse_z2

real
psc_p_pulse_z2(real x, real y, real z, real t)
{
  if (!psc.pulse_p_z2)
    return 0.;

  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_p_z2) { // default to Fortran
    return PSC_p_pulse_z2(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_p_z2, x, y, z, t);
}

//-----------------------------------------------------------------------
// psc_s_pulse_z2

real
psc_s_pulse_z2(real x, real y, real z, real t)
{
  // FIXME, create a fortran pulse instead of special casing
  if (!psc.pulse_s_z2) { // default to Fortran
    return PSC_s_pulse_z2(x, y, z, t);
  }
  return psc_pulse_field(psc.pulse_s_z2, x, y, z, t);
}

