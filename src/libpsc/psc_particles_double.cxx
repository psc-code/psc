
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_particles_inc.h"

// ======================================================================
// psc_mparticles: subclass "double"

template<>
mrc_obj_method MparticlesDouble::methods[] = {
  {}
};

#include "psc_particles_common.cxx"

