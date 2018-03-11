
#include "psc.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "double" particles are shifted that way.

#include "fields_item_moments_1st.hxx"

MAKE_POFI_OPS(PscMparticlesDouble, PscMfieldsC, double);
