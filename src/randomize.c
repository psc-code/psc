
#include "psc.h"
#include "util/profile.h"

// ======================================================================
// none randomize
//
// This doesn't actually randomize, it's used to turn randomize off.

static void
none_randomize()
{
}

struct psc_randomize_ops psc_randomize_ops_none = {
  .name      = "none",
  .randomize = none_randomize,
};
