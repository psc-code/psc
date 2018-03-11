
#include "psc.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single" particles are shifted that way.

#include "fields_item_moments_1st.hxx"

#include "psc_particles_single.h"
#include "psc_fields_c.h"

MAKE_POFI_OPS(PscMparticlesSingle, PscMfieldsC, single);
