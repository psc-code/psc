
#include "psc_particles_double.h"
#include "psc_particles_inc.h"

// ======================================================================
// psc_mparticles: subclass "double"

template<> const MparticlesBase::Convert MparticlesDouble::convert_to_ = {};
template<> const MparticlesBase::Convert MparticlesDouble::convert_from_ = {};

psc_mparticles_ops_<MparticlesDouble> psc_mparticles_double_ops;

