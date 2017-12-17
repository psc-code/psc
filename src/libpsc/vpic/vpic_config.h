
#ifndef VPIC_CONFIG_H
#define VPIC_CONFIG_H

#include "simulation.h"

#include "VpicGridBase.h"
#include "PscGridBase.h"

#include "VpicFieldArrayBase.h"
#include "PscFieldArrayBase.h"
#include "VpicFieldArrayLocalOps.h"
#include "PscFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "PscFieldArrayRemoteOps.h"
#include "VpicFieldArray.h"
#include "PscFieldArray.h"

#include "VpicParticleBc.h"

#include "VpicParticlesBase.h"
#include "PscParticlesBase.h"
#include "VpicParticlesOps.h"
#include "PscParticlesOps.h"

#include "VpicInterpolatorBase.h"
#include "PscInterpolatorBase.h"
#include "VpicInterpolator.h"
#include "PscInterpolator.h"

#include "VpicAccumulatorBase.h"
#include "PscAccumulatorBase.h"
#include "VpicAccumulator.h"
#include "PscAccumulator.h"

#include "VpicHydroArrayBase.h"
#include "PscHydroArrayBase.h"
#include "VpicHydroArray.h"
#include "PscHydroArray.h"

#include "VpicDiag.h"
#include "NoneDiag.h"

#include "VpicSimulationBase.h"
#include "PscSimulationBase.h"


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

typedef VpicParticleBcList ParticleBcList;

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
typedef VpicRng Rng;
typedef PscRngPool<Rng> RngPool;
#else
typedef VpicRng Rng;
typedef VpicRngPool<Rng> RngPool;
#endif

typedef VpicSimulation<Particles, ParticlesOps, RngPool, SimulationMixin, DiagMixin> Simulation;

#endif
