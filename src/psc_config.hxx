
#pragma once

template<typename DIM, typename Mparticles, typename Mfields>
struct PscConfigPushParticles2nd
{
  using PushParticles_t = PushParticles__<Config2nd<Mparticles, Mfields, DIM>>;
};

template<typename DIM, typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec
{
  using PushParticles_t = PushParticles1vb<Config1vbec<Mparticles, Mfields, DIM>>;
};

template<typename DIM, typename Mparticles, typename Mfields>
struct PscConfigPushParticlesCuda
{
};

template<typename Mparticles, typename Mfields>
struct PscConfigPushParticles1vbec<dim_xyz, Mparticles, Mfields>
{
  // need to use Config1vbecSplit when for dim_xyz
  using PushParticles_t = PushParticles1vb<Config1vbecSplit<Mparticles, Mfields, dim_xyz>>;
};

template<typename DIM, typename Mparticles, typename Mfields, template<typename...> class ConfigPushParticles>
struct PscConfig_
{
  using dim_t = DIM;
  using Mparticles_t = Mparticles;
  using Mfields_t = Mfields;
  using ConfigPushp = ConfigPushParticles<DIM, Mparticles, Mfields>;
  using PushParticles_t = typename ConfigPushp::PushParticles_t;
  using checks_order = typename PushParticles_t::checks_order;
  using Sort_t = SortCountsort2<Mparticles_t>;
  using Collision_t = Collision_<Mparticles_t, Mfields_t>;
  using PushFields_t = PushFields<Mfields_t>;
  using BndParticles_t = BndParticles_<Mparticles_t>;
  using Bnd_t = Bnd_<Mfields_t>;
  using BndFields_t = BndFieldsNone<Mfields_t>;
  using Inject_t = Inject_<Mparticles_t, MfieldsC, InjectFoil>; // FIXME, shouldn't always use MfieldsC
  using Heating_t = Heating__<Mparticles_t>;
  using Balance_t = Balance_<Mparticles_t, Mfields_t>;
  using Checks_t = Checks_<Mparticles_t, Mfields_t, checks_order>;
  using Marder_t = Marder_<Mparticles_t, Mfields_t>;
};

#ifdef USE_CUDA

template<typename DIM, typename Mparticles, typename Mfields>
struct PscConfig_<DIM, Mparticles, Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = DIM;
  using BS = typename Mparticles::BS;
  using Mparticles_t = Mparticles;
  using Mfields_t = Mfields;
  using PushParticles_t = PushParticlesCuda<CudaConfig1vbec3d<dim_t, BS>>;
  using Sort_t = SortCuda<BS>;
  using Collision_t = CollisionCuda<Mparticles>;
  using PushFields_t = PushFieldsCuda;
  using BndParticles_t = BndParticlesCuda<Mparticles, dim_t>;
  using Bnd_t = BndCuda3<Mfields_t>;
  using BndFields_t = BndFieldsNone<Mfields_t>;
  using Inject_t = InjectCuda<BS, dim_t, InjectFoil>;
  using Heating_t = HeatingCuda<BS>;
  using Balance_t = Balance_<MparticlesSingle, MfieldsSingle>;
  using Checks_t = ChecksCuda<Mparticles>;
  using Marder_t = MarderCuda<BS>;
};

template<typename Mparticles, typename Mfields>
struct PscConfig_<dim_xyz, Mparticles, Mfields, PscConfigPushParticlesCuda>
{
  using dim_t = dim_xyz;
  using BS = typename Mparticles::BS;
  using Mparticles_t = Mparticles;
  using Mfields_t = Mfields;
  using PushParticles_t = PushParticlesCuda<CudaConfig1vbec3dGmem<dim_t, BS>>;
  using Sort_t = SortCuda<BS>;
  using Collision_t = CollisionCuda<Mparticles>;
  using PushFields_t = PushFieldsCuda;
  using BndParticles_t = BndParticlesCuda<Mparticles, dim_t>;
  using Bnd_t = BndCuda3<Mfields>;
  using BndFields_t = BndFieldsNone<Mfields_t>;
  using Inject_t = InjectCuda<BS, dim_t, InjectFoil>;
  using Heating_t = HeatingCuda<BS>;
  using Balance_t = Balance_<MparticlesSingle, MfieldsSingle>;
  using Checks_t = ChecksCuda<Mparticles>;
  using Marder_t = MarderCuda<BS>;
};

#endif


template<typename dim>
using PscConfig2ndDouble = PscConfig_<dim, MparticlesDouble, MfieldsC, PscConfigPushParticles2nd>;

template<typename dim>
using PscConfig2ndSingle = PscConfig_<dim, MparticlesSingle, MfieldsSingle, PscConfigPushParticles2nd>;

template<typename dim>
using PscConfig1vbecSingle = PscConfig_<dim, MparticlesSingle, MfieldsSingle, PscConfigPushParticles1vbec>;

#ifdef USE_CUDA

template<typename dim>
struct PscConfig1vbecCuda : PscConfig_<dim, MparticlesCuda<BS144>, MfieldsCuda, PscConfigPushParticlesCuda>
{};

template<>
struct PscConfig1vbecCuda<dim_xyz> : PscConfig_<dim_xyz, MparticlesCuda<BS444>, MfieldsCuda, PscConfigPushParticlesCuda>
{};

#endif

