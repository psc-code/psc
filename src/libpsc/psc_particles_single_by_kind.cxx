
#include "bk_mparticles_iface.h"

#include "psc.h"
#include "psc_particles_single_by_kind.h"

#define PARTICLE_TYPE "single_by_kind"
#define PFX(x) psc_mparticles_single_by_kind_ ## x
#define psc_mparticles_sub psc_mparticles_single_by_kind

// ======================================================================
// psc_mparticles: subclass "single_by_kind"

// ----------------------------------------------------------------------
// psc_mparticles_ops

psc_mparticles_ops_<MparticlesSingleByKind> psc_mparticles_single_by_kind_ops;

