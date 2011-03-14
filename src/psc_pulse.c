
#include "psc.h"
#include "psc_glue.h"

// ----------------------------------------------------------------------
// psc_p_pulse_z1

real
psc_p_pulse_z1(real x, real y, real z, real t)
{
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

