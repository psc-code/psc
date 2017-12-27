
#ifndef VPIC_CONFIG_H
#define VPIC_CONFIG_H

//#define HAVE_VPIC

#ifndef HAVE_VPIC
//#include "util/profile/profile.h"
#undef TIC
#undef TOC
#define TIC
#define TOC(a, b)
#endif

#include "simulation.h"

#include "PscRng.h"

#include "PscGridBase.h"
#include "PscMaterial.h"

#include "PscFieldArrayBase.h"
#include "PscFieldArrayLocalOps.h"
#include "PscFieldArrayRemoteOps.h"
#include "PscFieldArray.h"

#include "PscParticleBc.h"

#include "PscParticlesBase.h"
#include "PscParticlesOps.h"

#include "PscInterpolatorBase.h"
#include "PscInterpolator.h"

#include "PscAccumulatorBase.h"
#include "PscAccumulator.h"

#include "PscHydroArrayBase.h"
#include "PscHydroArray.h"

#include "NoneDiag.h"

#include "PscSimulationBase.h"

#ifdef HAVE_VPIC

#include "VpicRng.h"

#include "VpicGridBase.h"

#include "VpicFieldArrayBase.h"
#include "VpicFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "VpicFieldArray.h"

#include "VpicParticleBc.h"

#include "VpicParticlesBase.h"
#include "VpicParticlesOps.h"

#include "VpicInterpolatorBase.h"
#include "VpicInterpolator.h"

#include "VpicAccumulatorBase.h"
#include "VpicAccumulator.h"

#include "VpicHydroArrayBase.h"
#include "VpicHydroArray.h"

#include "VpicDiag.h"

#include "VpicSimulationBase.h"
#endif


#if 1
typedef PscGridBase Grid;
#else
typedef VpicGridBase Grid;
#endif

#if 1
typedef PscMaterialList MaterialList;
typedef PscFieldArrayBase<Grid, MaterialList> FieldArrayBase;
typedef PscFieldArrayLocalOps<FieldArrayBase> FieldArrayLocalOps;
typedef PscFieldArrayRemoteOps<FieldArrayBase> FieldArrayRemoteOps;
typedef PscFieldArray<FieldArrayBase, FieldArrayLocalOps, FieldArrayRemoteOps> FieldArray;
#else
typedef VpicFieldArrayBase<Grid, VpicMaterialList> FieldArrayBase;
typedef VpicFieldArray<FieldArrayBase> FieldArray;
#endif

typedef PscInterpolatorBase<Grid> InterpolatorBase;
typedef PscInterpolator<InterpolatorBase, FieldArrayBase> Interpolator;

typedef PscAccumulatorBase<Grid> AccumulatorBase;
typedef PscAccumulator<AccumulatorBase, FieldArrayBase> Accumulator;

typedef PscHydroArrayBase<Grid> HydroArrayBase;
typedef PscHydroArray<HydroArrayBase> HydroArray;

#if 1
typedef PscParticleBcList ParticleBcList;
#else
typedef VpicParticleBcList ParticleBcList;
#endif

typedef PscParticlesBase<Grid, ParticleBcList> ParticlesBase;
typedef PscParticles<ParticlesBase, FieldArray, Interpolator, Accumulator, HydroArray> Particles;
#if 1
typedef PscParticlesOps<Particles> ParticlesOps;
#else
typedef VpicParticlesOps<Particles> ParticlesOps;
#endif

#if 1
typedef VpicDiagMixin<Particles> DiagMixin;
#else
typedef NoneDiagMixin<Particles> DiagMixin;
#endif

#if 1
typedef PscSimulationMixin<Particles> SimulationMixin;
#else
typedef VpicSimulationMixin<Particles> SimulationMixin;
#endif

#if 1
typedef PscRng Rng;
typedef PscRngPool<Rng> RngPool;
#else
typedef VpicRng Rng;
typedef VpicRngPool<Rng> RngPool;
#endif

typedef VpicSimulation<Particles, ParticlesOps, RngPool, SimulationMixin, DiagMixin> Simulation;

#endif

