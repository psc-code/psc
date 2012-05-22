
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

// ======================================================================
// psc_push_particles: subclass "1st"

struct psc_push_particles_ops psc_push_particles_1st_ops = {
  .name                  = "1st",
  .push_xz               = psc_push_particles_1st_push_xz,
  .push_a_yz             = psc_push_particles_1st_push_a_yz,
};

// ======================================================================
// psc_push_particles: subclass "1vb"

struct psc_push_particles_ops psc_push_particles_1vb_ops = {
  .name                  = "1vb",
  .push_a_yz             = psc_push_particles_1vb_push_a_yz,
};

// ======================================================================
// psc_push_particles: subclass "1sff"
//
// 1st order, self-force free

struct psc_push_particles_ops psc_push_particles_1sff_ops = {
  .name                  = "1sff",
  .push_xz               = psc_push_particles_1sff_push_xz,
};
