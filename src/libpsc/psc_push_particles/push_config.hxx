
#pragma once

#include "dim.hxx"
#include "fields.hxx"

#include "inc_defs.h"
#include "interpolate.hxx"
#include "inc_curr.c"
#include "pushp_current_esirkepov.hxx"
#include "inc_push.c"
#include "push_particles.hxx"

#define atomicAdd(addr, val) \
  do { *(addr) += (val); } while (0)

template<typename fields_t, typename dim_curr>
struct curr_cache_t : fields_t
{
  using real_t = typename fields_t::real_t;
  
  curr_cache_t(const fields_t& f)
    : fields_t(f.grid_, f.ib(), f.im(), f.n_comps(), f.data_)
  {}
  
  void add(int m, int i, int j, int k, real_t val)
  {
    Fields3d<fields_t, dim_curr> J(*this);
    real_t *addr = &J(JXI+m, i,j,k);
    atomicAdd(addr, val);
  }
};

template<typename _Mparticles, typename _MfieldsState,
	 typename _InterpolateEM,
	 typename _Dim, typename _Order,
	 class _Current>
struct _push_p_config
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using dim = _Dim;
  using Current_t = _Current;
  using InterpolateEM_t = _InterpolateEM;
  using AdvanceParticle_t = AdvanceParticle<typename Mparticles::real_t, _Dim>;
  using CurrentE_t = _Current;
};

template<typename _Mparticles, typename _MfieldsState,
	 typename _InterpolateEM,
	 typename _Dim, typename _Order>
using PushpConfigEsirkepov = _push_p_config<_Mparticles, _MfieldsState, _InterpolateEM,
					    _Dim, _Order,
					    Current<_Order, _Dim, Fields3d<typename _MfieldsState::fields_t>, _InterpolateEM>>;

template<typename _Mparticles, typename _MfieldsState,
	 typename _InterpolateEM,
	 typename _Dim, typename _Order,
	 template<typename, typename, typename> class _Current,
	 typename dim_curr = dim_xyz>
using PushpConfigVb = _push_p_config<_Mparticles, _MfieldsState, _InterpolateEM,
				     _Dim, _Order,
				     _Current<_Order, _Dim, curr_cache_t<typename _MfieldsState::fields_t, dim_curr>>>;

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

template<typename Mparticles, typename MfieldsState, typename dim>
using Config2nd = PushpConfigEsirkepov<Mparticles, MfieldsState,
				       InterpolateEM2nd<Fields3d<typename MfieldsState::fields_t>, dim>,
				       dim, opt_order_2nd>;

template<typename dim>
using Config2ndDouble = Config2nd<MparticlesDouble, MfieldsStateDouble, dim>;

template<typename dim>
using Config1stDouble = PushpConfigEsirkepov<MparticlesDouble, MfieldsStateDouble,
					     InterpolateEM1st<Fields3d<MfieldsStateDouble::fields_t>, dim>,
					     dim, opt_order_1st>;

template<typename Mparticles, typename Mfields, typename dim>
using Config1vbec = PushpConfigVb<Mparticles, Mfields,
				  InterpolateEM1vbec<Fields3d<typename Mfields::fields_t>, dim>,
				  dim, opt_order_1st,
				  Current1vbVar1>;

template<typename Mparticles, typename MfieldsState, typename dim>
using Config1vbecSplit = PushpConfigVb<Mparticles, MfieldsState,
				       InterpolateEM1vbec<Fields3d<typename MfieldsState::fields_t>, dim>,
				       dim, opt_order_1st,
				       Current1vbSplit>;

template<typename dim>
using Config1vbecDouble = Config1vbec<MparticlesDouble, MfieldsStateDouble, dim>;

template<typename dim>
using Config1vbecSingle = Config1vbec<MparticlesSingle, MfieldsStateSingle, dim>;

using Config1vbecSingleXZ = PushpConfigVb<MparticlesSingle, MfieldsStateSingle,
					  InterpolateEM1vbec<Fields3d<MfieldsStateSingle::fields_t, dim_xz>, dim_xyz>,
					  dim_xyz, opt_order_1st,
					  Current1vbSplit,
					  dim_xz>;
using Config1vbecSingle1 = PushpConfigVb<MparticlesSingle, MfieldsStateSingle,
					 InterpolateEM1vbec<Fields3d<MfieldsStateSingle::fields_t, dim_1>, dim_1>,
					 dim_1, opt_order_1st,
					 Current1vbVar1,
					 dim_1>;

