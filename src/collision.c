
#include "psc.h"
#include "util/profile.h"

// ======================================================================
// none collision
//
// This doesn't actually collision, it's used to turn collisions off.

static void
none_collision()
{
}

struct psc_collision_ops psc_collision_ops_none = {
  .name      = "none",
  .collision = none_collision,
};
