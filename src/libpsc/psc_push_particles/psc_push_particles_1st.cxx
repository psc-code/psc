
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

// ======================================================================
// psc_push_particles: subclass "1st"

struct psc_push_particles_ops_1st : psc_push_particles_ops {
  psc_push_particles_ops_1st() {
    name                  = "1st";
    push_mprts_xz         = psc_push_particles_push_mprts_1st_xz;
    push_mprts_yz         = psc_push_particles_push_mprts_1st_yz;
    particles_type        = "double";
  }
} psc_push_particles_1st_ops;

