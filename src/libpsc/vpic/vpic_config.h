
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

#include "PscFieldArrayLocalOps.h"
#include "PscFieldArrayRemoteOps.h"
#include "PscFieldArray.h"

#include "PscParticleBc.h"

#include "mparticles_vpic.hxx"
#include "PscParticlesBase.h"
#include "PscParticlesOps.h"

#include "mfields_interpolator_psc.hxx"
#include "PscInterpolator.h"

#include "PscAccumulator.h"
#include "mfields_accumulator_psc.hxx"

#include "mfields_hydro.hxx"

#include "NoneDiag.h"

struct VpicConfigPsc
{
  using Grid = PscGridBase;

  using ParticleBcList = PscParticleBcList;
  using Particles = PscParticlesBase<Grid, ParticleBcList>;
  using Mparticles = MparticlesVpic_<Particles>;

  using MaterialList = PscMaterialList;
  using MfieldsState = MfieldsStatePsc<Grid, MaterialList>;

  using FieldArrayLocalOps = PscFieldArrayLocalOps<MfieldsState>;
  using FieldArrayRemoteOps = PscFieldArrayRemoteOps<MfieldsState>;
  using AccumulateOps =
    PscAccumulateOps<MfieldsState, FieldArrayLocalOps, FieldArrayRemoteOps>;

  using MfieldsInterpolator = MfieldsInterpolatorPsc<Grid>;
  using InterpolatorOps = PscInterpolatorOps<MfieldsInterpolator, MfieldsState>;

  using MfieldsAccumulator = MfieldsAccumulatorPsc<Grid>;
  using AccumulatorOps = PscAccumulatorOps<MfieldsAccumulator, MfieldsState>;

  using MfieldsHydro = MfieldsHydroPsc<Grid>;

  using ParticlesOps =
    PscParticlesOps<Mparticles, MfieldsState, MfieldsInterpolator,
                    MfieldsAccumulator, MfieldsHydro>;

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
