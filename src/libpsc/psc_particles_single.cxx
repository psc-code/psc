
#include "psc.h"
#include "psc_particles_single.h"
#include "psc_particles_inc.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_mparticles: subclass "single"

// ----------------------------------------------------------------------
// conversion to/from "double"

template<> const MparticlesBase::Convert MparticlesSingle::convert_to_ = {
  { std::type_index(typeid(MparticlesDouble)), psc_mparticles_copy_to<MparticlesSingle, MparticlesDouble> },
};

template<> const MparticlesBase::Convert MparticlesSingle::convert_from_ = {
  { std::type_index(typeid(MparticlesDouble)), psc_mparticles_copy_from<MparticlesSingle, MparticlesDouble> },
};

psc_mparticles_ops_<MparticlesSingle> psc_mparticles_single_ops;
