
#include "psc_randomize_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_randomize_none_run

static void
psc_randomize_none_run(struct psc_randomize *randomize,
		       struct psc_mparticles *mprts_base)
{
}

// ======================================================================
// psc_randomize: subclass "none"

struct psc_randomize_ops_none : psc_randomize_ops {
  psc_randomize_ops_none() {
    name                  = "none";
    run                   = psc_randomize_none_run;
  }
} psc_randomize_none_ops;
