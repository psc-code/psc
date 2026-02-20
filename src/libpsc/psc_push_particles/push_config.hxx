
#pragma once

#include "dim.hxx"
#include "fields.hxx"

#include "inc_defs.h"
#include "interpolate.hxx"
#include "inc_curr.cxx"
#include "inc_push.cxx"
#include "push_particles.hxx"
#include "push_particles_esirkepov.hxx"
#include "push_particles_1vb.hxx"

#include <psc/gtensor.h>

template <typename fields_t>
class curr_cache_t
{
public:
  using real_t = typename fields_t::value_type;
  using value_type = typename fields_t::value_type;
  using storage_type = typename fields_t::Storage;

  curr_cache_t(fields_t& f) : storage_(f.storage()), ib_(f.ib()) {}

  void add(int m, int i, int j, int k, real_t val)
  {
    storage_(i - ib_[0], j - ib_[1], k - ib_[2], JXI + m) += val;
  }

private:
  storage_type storage_;
  Int3 ib_;
};

template <typename _Mparticles, typename _MfieldsState,
          template <typename, typename> class _InterpolateEM, typename _Dim,
          typename _Order>
struct PushpConfigEsirkepov
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Dim = _Dim;
  using Order = _Order;
  using InterpolateEM_t =
    _InterpolateEM<Fields3d<typename MfieldsState::fields_view_t::Storage>,
                   Dim>;
  using AdvanceParticle_t = AdvanceParticle<typename Mparticles::real_t, Dim>;
};

template <typename _Mparticles, typename _MfieldsState, typename _InterpolateEM,
          typename _Dim, typename _Order,
          template <typename, typename, typename> class _Current>
struct PushpConfigVb
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Dim = _Dim;
  using InterpolateEM_t = _InterpolateEM;
  using Current_t =
    _Current<_Order, _Dim, curr_cache_t<typename _MfieldsState::fields_view_t>>;
  using AdvanceParticle_t = AdvanceParticle<typename Mparticles::real_t, Dim>;
};

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

template <typename Mparticles, typename MfieldsState, typename dim>
using Config2nd = PushpConfigEsirkepov<Mparticles, MfieldsState,
                                       InterpolateEM2nd, dim, opt_order_2nd>;

template <typename dim>
using Config2ndDouble = Config2nd<MparticlesDouble, MfieldsStateDouble, dim>;

template <typename dim>
using Config1stDouble =
  PushpConfigEsirkepov<MparticlesDouble, MfieldsStateDouble, InterpolateEM1st,
                       dim, opt_order_1st>;

template <typename Mparticles, typename Mfields, typename dim>
using Config1vbec = PushpConfigVb<
  Mparticles, Mfields,
  InterpolateEM1vbec<Fields3d<typename Mfields::fields_view_t::Storage>, dim>,
  dim, opt_order_1st, Current1vbVar1>;

template <typename Mparticles, typename MfieldsState, typename dim>
using Config1vbecSplit =
  PushpConfigVb<Mparticles, MfieldsState,
                InterpolateEM1vbec<
                  Fields3d<typename MfieldsState::fields_view_t::Storage>, dim>,
                dim, opt_order_1st, Current1vbSplit>;

template <typename dim>
using Config1vbecDouble =
  Config1vbec<MparticlesDouble, MfieldsStateDouble, dim>;

template <typename dim>
using Config1vbecSingle =
  Config1vbec<MparticlesSingle, MfieldsStateSingle, dim>;

using Config1vbecSingleXZ = PushpConfigVb<
  MparticlesSingle, MfieldsStateSingle,
  InterpolateEM1vbec<
    Fields3d<MfieldsStateSingle::fields_view_t::Storage, dim_xz>, dim_xyz>,
  dim_xyz, opt_order_1st, Current1vbSplit>;
using Config1vbecSingle1 = PushpConfigVb<
  MparticlesSingle, MfieldsStateSingle,
  InterpolateEM1vbec<
    Fields3d<MfieldsStateSingle::fields_view_t::Storage, dim_1>, dim_1>,
  dim_1, opt_order_1st, Current1vbVar1>;
