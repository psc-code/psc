
#pragma once

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"

#include "../libpsc/psc_sort/psc_sort_impl.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"
#include "../libpsc/psc_push_particles/push_part_common.c"
#include "../libpsc/psc_push_particles/1vb/psc_push_particles_1vb.h"
#include "../libpsc/psc_push_particles/1vb.c"
#include "psc_push_fields_impl.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "../libpsc/psc_bnd_fields/psc_bnd_fields_impl.hxx"
#include "bnd_particles_impl.hxx"
#include "../libpsc/psc_balance/psc_balance_impl.hxx"
#include "../libpsc/psc_push_fields/marder_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/push_particles_cuda_impl.hxx"
#include "../libpsc/cuda/push_fields_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_2_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_3_impl.hxx"
#include "../libpsc/cuda/bnd_particles_cuda_impl.hxx"
#include "../libpsc/cuda/sort_cuda_impl.hxx"
#include "../libpsc/cuda/collision_cuda_impl.hxx"
#include "../libpsc/cuda/checks_cuda_impl.hxx"
#include "../libpsc/cuda/marder_cuda_impl.hxx"
#endif

struct SimulationNone
{
  using Species = void;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticles2nd
{
  using PushParticles_t = PushParticles__<Config2nd<Mparticles, MfieldsState, DIM>>;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticles1vbec
{
  using PushParticles_t = PushParticles1vb<Config1vbec<Mparticles, MfieldsState, DIM>>;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticlesCuda
{
};

// need to use Config1vbecSplit when for dim_xyz, dim_xz

template<typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec<dim_xyz, Mparticles, Mfields>
{
  using PushParticles_t = PushParticles1vb<Config1vbecSplit<Mparticles, Mfields, dim_xyz>>;
};

template<typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec<dim_xz, Mparticles, Mfields>
{
  using PushParticles_t = PushParticles1vb<Config1vbecSplit<Mparticles, Mfields, dim_xz>>;
};

template<typename DIM, typename Mparticles, typename _MfieldsState,
	 typename _Mfields,
	 template<typename...> class ConfigPushParticles,
	 typename _Simulation = SimulationNone>
struct PscConfig_
{
  using dim_t = DIM;
  using Mparticles_t = Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using ConfigPushp = ConfigPushParticles<DIM, Mparticles, MfieldsState>;
  using PushParticles_t = typename ConfigPushp::PushParticles_t;
  using checks_order = typename PushParticles_t::checks_order;
  using Sort_t = SortCountsort2<Mparticles_t>;
  using Collision_t = Collision_<Mparticles_t, MfieldsState, Mfields>;
  using PushFields_t = PushFields<MfieldsState>;
  using BndParticles_t = BndParticles_<Mparticles_t>;
  using Bnd_t = Bnd_<MfieldsState>;
  using BndFields_t = BndFieldsNone<MfieldsState>;
  using Balance_t = Balance_<Mparticles_t, MfieldsState, Mfields>;
  using Checks_t = Checks_<Mparticles_t, MfieldsState, Mfields, checks_order>;
  using Marder_t = Marder_<Mparticles_t, MfieldsState, Mfields>;
  using Simulation = _Simulation;
};

#ifdef USE_CUDA

template<typename DIM, typename Mparticles, typename _MfieldsState, typename _Mfields>
struct PscConfig_<DIM, Mparticles, _MfieldsState, _Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = DIM;
  using BS = typename Mparticles::BS;
  using Mparticles_t = Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using PushParticles_t = PushParticlesCuda<CudaConfig1vbec3d<dim_t, BS>>;
  using Sort_t = SortCuda<BS>;
  using Collision_t = CollisionCuda<Mparticles>;
  using PushFields_t = PushFieldsCuda;
  using BndParticles_t = BndParticlesCuda<Mparticles, dim_t>;
  using Bnd_t = BndCuda3<MfieldsState>;
  using BndFields_t = BndFieldsNone<MfieldsState>;
  using Balance_t = Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
  using Checks_t = ChecksCuda<Mparticles>;
  using Marder_t = MarderCuda<BS>;
};

template<typename Mparticles, typename _MfieldsState, typename _Mfields>
struct PscConfig_<dim_xyz, Mparticles, _MfieldsState, _Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = dim_xyz;
  using BS = typename Mparticles::BS;
  using Mparticles_t = Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using PushParticles_t = PushParticlesCuda<CudaConfig1vbec3dGmem<dim_t, BS>>;
  using Sort_t = SortCuda<BS>;
  using Collision_t = CollisionCuda<Mparticles>;
  using PushFields_t = PushFieldsCuda;
  using BndParticles_t = BndParticlesCuda<Mparticles, dim_t>;
  using Bnd_t = BndCuda3<MfieldsState>;
  using BndFields_t = BndFieldsNone<MfieldsState>;
  using Balance_t = Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
  using Checks_t = ChecksCuda<Mparticles>;
  using Marder_t = MarderCuda<BS>;
};

#endif


template<typename dim>
using PscConfig2ndDouble = PscConfig_<dim, MparticlesDouble, MfieldsStateDouble, MfieldsC,
				      PscConfigPushParticles2nd>;

template<typename dim>
using PscConfig2ndSingle = PscConfig_<dim, MparticlesSingle, MfieldsStateSingle, MfieldsSingle,
				      PscConfigPushParticles2nd>;

template<typename dim>
using PscConfig1vbecSingle = PscConfig_<dim, MparticlesSingle, MfieldsStateSingle, MfieldsSingle,
					PscConfigPushParticles1vbec>;

template<typename dim>
using PscConfig1vbecDouble = PscConfig_<dim, MparticlesDouble, MfieldsStateDouble, MfieldsC,
					PscConfigPushParticles1vbec>;

#ifdef USE_CUDA

template<typename dim>
struct PscConfig1vbecCuda : PscConfig_<dim, MparticlesCuda<BS144>, MfieldsStateCuda, MfieldsCuda, PscConfigPushParticlesCuda>
{};

template<>
struct PscConfig1vbecCuda<dim_xyz> : PscConfig_<dim_xyz, MparticlesCuda<BS444>, MfieldsStateCuda, MfieldsCuda, PscConfigPushParticlesCuda>
{};

#endif

#ifdef USE_VPIC

#include "../libpsc/vpic/sort_vpic.hxx"
#include "../libpsc/vpic/collision_vpic.hxx"
#include "../libpsc/vpic/push_particles_vpic.hxx"
#include "../libpsc/vpic/push_fields_vpic.hxx"
#include "../libpsc/vpic/bnd_vpic.hxx"
#include "../libpsc/vpic/bnd_fields_vpic.hxx"
#include "../libpsc/vpic/bnd_particles_vpic.hxx"
#include "../libpsc/vpic/marder_vpic.hxx"
#include "../libpsc/vpic/checks_vpic.hxx"

struct PscConfigVpic
{
  using Mparticles_t = MparticlesVpic;
  using MfieldsState = MfieldsState;
  using Balance_t = Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
  using Sort_t = SortVpic;
  using Collision_t = PscCollisionVpic;
  using PushParticles_t = PushParticlesVpic;
  using PushFields_t = PushFieldsVpic;
  using Bnd_t = BndVpic;
  using BndFields_t = BndFieldsVpic;
  using BndParticles_t = BndParticlesVpic;
  using Checks_t = ChecksVpic;
  using Marder_t = MarderVpic;

#if 0
  using DiagMixin = VpicDiagMixin<MparticlesVpic, MfieldsState, MfieldsInterpolator, MfieldsHydro,
				  DiagOps, ParticlesOps, HydroArrayOps>;
#else
  using DiagMixin = NoneDiagMixin<MparticlesVpic, MfieldsState, MfieldsInterpolator, MfieldsHydro>;
#endif
};

#endif
