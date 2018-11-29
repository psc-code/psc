
#ifndef VPIC_CONFIG_H
#define VPIC_CONFIG_H

#define DO_VPIC 1

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

#endif

struct VpicConfig_
{
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
  
#ifdef DO_VPIC
  using PushFieldsOps = VpicPushFieldsOps<MfieldsState>;
  using AccumulateOps = VpicAccumulateOps<MfieldsState>;
  using CleanDivOps = VpicCleanDivOps<MfieldsState>;
#else
  using FieldArrayLocalOps = PscFieldArrayLocalOps<MfieldsState>;
  using FieldArrayRemoteOps = PscFieldArrayRemoteOps<MfieldsState>;
  using PushFieldsOps = PscPushFieldsOps<MfieldsState, FieldArrayLocalOps, FieldArrayRemoteOps>;
  using AccumulateOps = PscAccumulateOps<MfieldsState, FieldArrayLocalOps, FieldArrayRemoteOps>;
  using CleanDivOps = PscCleanDivOps<MfieldsState, FieldArrayLocalOps, FieldArrayRemoteOps>;
#endif

#ifdef DO_VPIC
  using MfieldsInterpolator = MfieldsInterpolatorVpic;
#else
  using MfieldsInterpolator = MfieldsInterpolatorPsc<Grid>;
#endif
  using InterpolatorOps = PscInterpolatorOps<MfieldsInterpolator, MfieldsState>;
  
#ifdef DO_VPIC
  using MfieldsAccumulator = MfieldsAccumulatorVpic;
#else
  using MfieldsAccumulator = MfieldsAccumulatorPsc<Grid>;
#endif
  using AccumulatorOps = PscAccumulatorOps<MfieldsAccumulator, MfieldsState>;
  
#ifdef DO_VPIC
  using MfieldsHydro = MfieldsHydroVpic;
  using HydroArrayOps = VpicHydroArrayOps<MfieldsHydro>;
#else
  using MfieldsHydro = MfieldsHydroPsc<Grid>;
  using HydroArrayOps = PscHydroArrayOps<MfieldsHydro>;
#endif

#ifdef DO_VPIC
  using ParticleBcList = VpicParticleBcList;
  using Particles = VpicParticlesBase<Grid, ParticleBcList>;
#else
  using ParticleBcList = PscParticleBcList;
  using Particles = PscParticlesBase<Grid, ParticleBcList>;
#endif
  using Mparticles = MparticlesVpic_<Particles>;

#if 0//def DO_VPIC
  using ParticlesOps = VpicParticlesOps<Particles, MfieldsState, Interpolator, MfieldsAccumulator, MfieldsHydro>;
#else
  using ParticlesOps = PscParticlesOps<Mparticles, MfieldsState, MfieldsInterpolator, MfieldsAccumulator, MfieldsHydro>;
#endif

/* #ifdef DO_VPIC */
/*   using DiagOps = VpicDiagOps<MfieldsState>; */
/* #else */
/*   using DiagOps = PscDiagOps<MfieldsState>; */
/* #endif */
};

#if 1
typedef PscRng Rng;
typedef PscRngPool<Rng> RngPool;
#else
typedef VpicRng Rng;
typedef VpicRngPool<Rng> RngPool;
#endif

#endif

