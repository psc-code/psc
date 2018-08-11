
#ifndef VPIC_CONFIG_H
#define VPIC_CONFIG_H

#define HAVE_VPIC

#ifdef HAVE_VPIC
#include "util/profile/profile.h"
#else
//#include "util/profile/profile.h"
//#undef TIC
//#undef TOC
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
#include "PscPushFieldsOps.hxx"

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

#define particle_t __vpic_particle_t

#include "VpicRng.h"

#include "VpicGridBase.h"
#include "VpicMaterial.h"

#include "VpicFieldArrayBase.h"
#include "VpicFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "VpicFieldArray.h"
#include "VpicPushFieldsOps.hxx"

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

#undef particle_t

#endif


#if 1
typedef PscGridBase Grid;
#else
typedef VpicGridBase Grid;
#endif

#if 1
using MaterialList = PscMaterialList;
using FieldArray = PscFieldArrayBase<Grid, MaterialList>;
using FieldArrayLocalOps = PscFieldArrayLocalOps<FieldArray>;
using FieldArrayRemoteOps = PscFieldArrayRemoteOps<FieldArray>;
using PushFieldsOps = PscPushFieldsOps<FieldArray, FieldArrayLocalOps, FieldArrayRemoteOps>;
using DiagOps = PscDiagOps<FieldArray>;
using AccumulateOps = PscAccumulateOps<FieldArray, FieldArrayLocalOps, FieldArrayRemoteOps>;
using CleanDivOps = PscCleanDivOps<FieldArray, FieldArrayLocalOps, FieldArrayRemoteOps>;
#else
using FieldArray = VpicFieldArrayBase<Grid, VpicMaterialList>;
using PushFieldsOps = VpicPushFieldsOps<FieldArray>;
using DiagOps = VpicDiagOps<FieldArray>;
using AccumulateOps = VpicAccumulateOps<FieldArray>;
using CleanDivOps = VpicCleanDivOps<FieldArray>;
#endif

using Interpolator = PscInterpolatorBase<Grid>;
using InterpolatorOps = PscInterpolatorOps<Interpolator, FieldArray>;

using Accumulator = PscAccumulatorBase<Grid>;
using AccumulatorOps = PscAccumulatorOps<Accumulator, FieldArray>;

using HydroArray = PscHydroArrayBase<Grid>;
using HydroArrayOps = PscHydroArrayOps<HydroArray>;

#if 1
using ParticleBcList = PscParticleBcList;
#else
using ParticleBcList = VpicParticleBcList;
#endif

typedef PscParticlesBase<Grid, ParticleBcList> ParticlesBase;
typedef PscParticles<ParticlesBase, FieldArray, Interpolator, Accumulator, HydroArray> Particles;
#if 1
typedef PscParticlesOps<Particles> ParticlesOps;
#else
typedef VpicParticlesOps<Particles> ParticlesOps;
#endif

#if 1
typedef VpicDiagMixin<Particles, DiagOps, HydroArrayOps> DiagMixin;
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

typedef VpicSimulation<Particles, RngPool, SimulationMixin, DiagMixin> Simulation;

#endif

