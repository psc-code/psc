
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

// ======================================================================
// psc_push_particles: subclass "1st"

struct psc_push_particles_ops psc_push_particles_1st_ops = {
  .name                  = "1st",
  .push_a_xz             = psc_push_particles_1st_push_a_xz,
  .push_mprts_yz         = psc_push_particles_1st_push_mprts_yz,
  .particles_type        = "c",
  .fields_type           = "c",
};

// ======================================================================
// psc_push_particles: subclass "1sff"
//
// 1st order, self-force free

struct psc_push_particles_ops psc_push_particles_1sff_ops = {
  .name                  = "1sff",
  .push_a_xz             = psc_push_particles_1sff_push_a_xz,
  .particles_type        = "c",
  .fields_type           = "c",
};
