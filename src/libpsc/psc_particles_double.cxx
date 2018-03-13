
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_particles_inc.h"

// ======================================================================
// psc_mparticles: subclass "double"

template<> const MparticlesBase::Map MparticlesDouble::convert_to_ = {};
template<> const MparticlesBase::Map MparticlesDouble::convert_from_ = {};

#include "psc_particles_common.cxx"

