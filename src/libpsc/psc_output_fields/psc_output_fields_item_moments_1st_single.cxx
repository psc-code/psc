
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single" particles are shifted that way.

#include "psc_output_fields_item_moments_1st.cxx"

MAKE_POFI_OPS(single);
