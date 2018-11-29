
#ifndef VPIC_CONFIG_H
#define VPIC_CONFIG_H

//#define DO_VPIC 1

#ifdef USE_VPIC
#include "util/profile/profile.h"
#else
//#include "util/profile/profile.h"
//#undef TIC
//#undef TOC
#define TIC
#define TOC(a, b)
#endif

#include "PscRng.h"

#include "PscGridBase.h"
#include "PscMaterial.h"

#include "mfields_state_psc.hxx"
#include "PscFieldArrayLocalOps.h"
#include "PscFieldArrayRemoteOps.h"
#include "PscFieldArray.h"
#include "PscPushFieldsOps.hxx"

#include "PscParticleBc.h"

#include "mparticles_vpic.hxx"
#include "PscParticlesBase.h"
#include "PscParticlesOps.h"

#include "mfields_interpolator_psc.hxx"
#include "PscInterpolator.h"

#include "PscAccumulator.h"
#include "mfields_accumulator_psc.hxx"

#include "PscHydroArray.h"
#include "mfields_hydro.hxx"

#include "NoneDiag.h"

#ifdef USE_VPIC

#include "VpicRng.h"

#include "VpicGridBase.h"
#include "VpicMaterial.h"

#include "VpicFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "VpicFieldArray.h"
#include "VpicPushFieldsOps.hxx"

#include "VpicParticleBc.h"

#include "VpicParticlesBase.h"
#include "VpicParticlesOps.h"

#include "VpicInterpolator.h"

#include "VpicAccumulator.h"

#include "VpicHydroArray.h"

#include "VpicDiag.h"

#include "mfields_state_vpic.hxx"
#include "mfields_hydro_vpic.hxx"
#include "mfields_interpolator_vpic.hxx"
#include "mfields_accumulator_vpic.hxx"

#undef particle_t

#endif

#ifdef DO_VPIC
using Grid = VpicGridBase;
#else
using Grid = PscGridBase;
#endif

#ifdef DO_VPIC
using MaterialList = VpicMaterialList;
using MfieldsState = MfieldsStateVpic;
#else
using MaterialList = PscMaterialList;
using MfieldsState = MfieldsStatePsc<Grid, MaterialList>;
#endif

#if 1
typedef PscRng Rng;
typedef PscRngPool<Rng> RngPool;
#else
typedef VpicRng Rng;
typedef VpicRngPool<Rng> RngPool;
#endif

#endif

