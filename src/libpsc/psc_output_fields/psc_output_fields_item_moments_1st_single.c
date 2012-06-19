
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single" particles are shifted that way.

#include "psc_output_fields_item_moments_1st.c"

MAKE_POFI_OPS(n_1st, single);
MAKE_POFI_OPS(v_1st, single);
MAKE_POFI_OPS(p_1st, single);
MAKE_POFI_OPS(vv_1st, single);
