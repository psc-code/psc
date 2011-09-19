
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

// ======================================================================
// psc_push_particles: subclass "1st"

struct psc_push_particles_ops psc_push_particles_1st_ops = {
  .name                  = "1st",
  .push_xz               = psc_push_particles_1st_push_xz,
};
