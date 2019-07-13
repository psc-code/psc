
#pragma once

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"

#include "../libpsc/psc_sort/psc_sort_impl.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"
#include "../libpsc/psc_push_particles/1vb/psc_push_particles_1vb.h"
#include "psc_push_fields_impl.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "../libpsc/psc_bnd_fields/psc_bnd_fields_impl.hxx"
#include "bnd_particles_impl.hxx"
#include "../libpsc/psc_balance/psc_balance_impl.hxx"
#include "../libpsc/psc_push_fields/marder_impl.hxx"
#include "../libpsc/psc_output_particles/output_particles_hdf5_impl.hxx"
#include "../libpsc/psc_output_particles/output_particles_none_impl.hxx"

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

template<typename Mparticles>
using OutputParticlesDefault = OutputParticlesHdf5<Mparticles>;

struct SimulationNone
{
  using Species = void;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticles2nd
{
  using PushParticles_t = PushParticlesEsirkepov<Config2nd<Mparticles, MfieldsState, DIM>>;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticles1vbec
{
  using PushParticles_t = PushParticlesVb<Config1vbec<Mparticles, MfieldsState, DIM>>;
};

template<typename DIM, typename Mparticles, typename MfieldsState>
struct PscConfigPushParticlesCuda
{
};

// need to use Config1vbecSplit when for dim_xyz, dim_xz

template<typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec<dim_xyz, Mparticles, Mfields>
{
  using PushParticles_t = PushParticlesVb<Config1vbecSplit<Mparticles, Mfields, dim_xyz>>;
};

template<typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec<dim_xz, Mparticles, Mfields>
{
  using PushParticles_t = PushParticlesVb<Config1vbecSplit<Mparticles, Mfields, dim_xz>>;
};

template<typename DIM, typename _Mparticles, typename _MfieldsState,
	 typename _Mfields,
	 template<typename...> class ConfigPushParticles,
	 typename _Simulation = SimulationNone>
struct PscConfig_
{
  using dim_t = DIM;
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using ConfigPushp = ConfigPushParticles<DIM, Mparticles, MfieldsState>;
  using PushParticles_t = typename ConfigPushp::PushParticles_t;
  using checks_order = typename PushParticles_t::checks_order;
  using Sort_t = SortCountsort2<Mparticles>;
  using Collision_t = Collision_<Mparticles, MfieldsState, Mfields>;
  using PushFields_t = PushFields<MfieldsState>;
  using BndParticles_t = BndParticles_<Mparticles>;
  using Bnd_t = Bnd_<MfieldsState>;
  using BndFields_t = BndFieldsNone<MfieldsState>;
  using Balance_t = Balance_<Mparticles, MfieldsState, Mfields>;
  using Checks_t = Checks_<Mparticles, MfieldsState, Mfields, checks_order>;
  using Marder_t = Marder_<Mparticles, MfieldsState, Mfields>;
  using Simulation = _Simulation;
  using OutputParticles = OutputParticlesDefault<Mparticles>;
};

#ifdef USE_CUDA

template<typename DIM, typename _Mparticles, typename _MfieldsState, typename _Mfields>
struct PscConfig_<DIM, _Mparticles, _MfieldsState, _Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = DIM;
  using BS = typename _Mparticles::BS;
  using Mparticles = _Mparticles;
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
  using OutputParticles = OutputParticlesDefault<MparticlesSingle>;
};

template<typename _Mparticles, typename _MfieldsState, typename _Mfields>
struct PscConfig_<dim_xyz, _Mparticles, _MfieldsState, _Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = dim_xyz;
  using BS = typename _Mparticles::BS;
  using Mparticles = _Mparticles;
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
  using OutputParticles = OutputParticlesNone;
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

#include "../libpsc/vpic/sort_vpic.hxx"
#include "../libpsc/vpic/collision_vpic.hxx"
#include "../libpsc/vpic/push_particles_vpic.hxx"
#include "../libpsc/vpic/push_fields_vpic.hxx"
#include "../libpsc/vpic/bnd_vpic.hxx"
#include "../libpsc/vpic/bnd_fields_vpic.hxx"
#include "../libpsc/vpic/bnd_particles_vpic.hxx"
#include "../libpsc/vpic/marder_vpic.hxx"
#include "../libpsc/vpic/checks_vpic.hxx"

#ifdef USE_VPIC
struct PscConfigVpicWrap
{
  using VpicConfig = VpicConfigWrap;
  
  using MfieldsState = typename VpicConfig::MfieldsState;
  using Mparticles = typename VpicConfig::Mparticles;
  using MfieldsHydro = typename VpicConfig::MfieldsHydro;

  using Balance_t = Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
  using Sort_t = SortVpicWrap<Mparticles>;
  using Collision_t = PscCollisionVpic;
  using PushParticles_t = PushParticlesVpic<Mparticles, MfieldsState,
					    typename VpicConfig::ParticlesOps,
					    typename VpicConfig::AccumulatorOps,
					    typename VpicConfig::AccumulateOps,
					    typename VpicConfig::InterpolatorOps>;
  using PushFields_t = PushFieldsVpicWrap<MfieldsState>;
  using Bnd_t = BndVpic<MfieldsState>;
  using BndFields_t = BndFieldsVpic<MfieldsState>;
  using BndParticles_t = BndParticlesVpic<Mparticles>;
  using Checks_t = ChecksVpic<Mparticles, MfieldsState>;
  using Marder_t = MarderVpicWrap<Mparticles, MfieldsState>;
  using OutputParticles = OutputParticlesHdf5<MparticlesSingle>;
  using OutputHydro = OutputHydroVpicWrap<Mparticles, MfieldsHydro, typename VpicConfig::MfieldsInterpolator>;
  using dim_t = dim_xyz;

#if 0
  using DiagMixin = VpicDiagMixin<Mparticles, MfieldsState, MfieldsInterpolator, MfieldsHydro,
				  DiagOps, ParticlesOps, HydroArrayOps>;
#else
  using DiagMixin = NoneDiagMixin<Mparticles, MfieldsState,
				  typename VpicConfig::MfieldsInterpolator,
				  typename VpicConfig::MfieldsHydro>;
#endif
};
#endif

struct PscConfigVpicPsc
{
  using VpicConfig = VpicConfigPsc;

  using MfieldsState = typename VpicConfig::MfieldsState;
  using Mparticles = typename VpicConfig::Mparticles;
  using MfieldsHydro = typename VpicConfig::MfieldsHydro;
  
  using Balance_t = Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>;
  using Sort_t = SortVpic<Mparticles>;
  using Collision_t = PscCollisionVpic;
  using PushParticles_t = PushParticlesVpic<Mparticles, MfieldsState,
					    typename VpicConfig::ParticlesOps,
					    typename VpicConfig::AccumulatorOps,
					    typename VpicConfig::AccumulateOps,
					    typename VpicConfig::InterpolatorOps>;
  using PushFields_t = PushFieldsVpic<MfieldsState>;
  using Bnd_t = BndVpic<MfieldsState>;
  using BndFields_t = BndFieldsVpic<MfieldsState>;
  using BndParticles_t = BndParticlesVpic<Mparticles>;
  using Checks_t = ChecksVpic<Mparticles, MfieldsState>;
  using Marder_t = MarderVpic<Mparticles, MfieldsState>;
  using OutputParticles = OutputParticlesHdf5<MparticlesSingle>;
  using OutputHydro = OutputHydroVpic<Mparticles, MfieldsHydro, typename VpicConfig::MfieldsInterpolator>;
  using dim_t = dim_xyz;

#if 0
  using DiagMixin = VpicDiagMixin<Mparticles, MfieldsState, MfieldsInterpolator, MfieldsHydro,
				  DiagOps, ParticlesOps, HydroArrayOps>;
#else
  using DiagMixin = NoneDiagMixin<Mparticles, MfieldsState,
				  typename VpicConfig::MfieldsInterpolator,
				  typename VpicConfig::MfieldsHydro>;
#endif
};

