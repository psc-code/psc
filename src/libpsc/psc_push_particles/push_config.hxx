
#pragma once

#include "dim.hxx"
#include "fields.hxx"

#include "inc_defs.h"
#include "pushp_cache_fields.hxx"
#include "interpolate.hxx"
#include "inc_curr.c"
#include "pushp_current_esirkepov.hxx"
#include "inc_push.c"

#define atomicAdd(addr, val) \
  do { *(addr) += (val); } while (0)

template<typename fields_t, typename dim_curr>
struct curr_cache_t : fields_t
{
  using real_t = typename fields_t::real_t;
  
  curr_cache_t(const fields_t& f)
    : fields_t(f.ib(), f.im(), f.n_comps(), f.data_)
  {}
  
  void add(int m, int i, int j, int k, real_t val)
  {
    Fields3d<fields_t, dim_curr> J(*this);
    real_t *addr = &J(JXI+m, i,j,k);
    atomicAdd(addr, val);
  }
};

template<typename MP, typename MF,
	 typename IPEM,
	 typename D, typename IP, typename O,
	 template<typename, typename> class CURRENT,
	 template<typename, typename> class CF = CacheFieldsNone,
	 typename dim_curr = dim_xyz>
struct push_p_config
{
  using Mparticles = MP;
  using Mfields = MF;
  using dim = D;
  using CacheFields = CF<typename MF::fields_t, D>;
  using curr_cache_t = curr_cache_t<typename MF::fields_t, dim_curr>;
  using Current_t = CURRENT<curr_cache_t, D>;

  using InterpolateEM_t = IPEM;
  using CurrentE_t = Current<O, D, InterpolateEM_t, Fields3d<typename MF::fields_t>>;
};

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

template<typename dim>
using Config2nd = push_p_config<MparticlesDouble, MfieldsC,
				InterpolateEM<Fields3d<MfieldsC::fields_t>, opt_ip_2nd, dim>,
				dim, opt_ip_2nd, opt_order_2nd,
				CurrentNone>;

using Config2ndDoubleYZ = push_p_config<MparticlesDouble, MfieldsC,
					InterpolateEM<Fields3d<MfieldsC::fields_t>, opt_ip_2nd, dim_yz>,
					dim_yz, opt_ip_2nd, opt_order_2nd,
					CurrentNone, CacheFields>;

template<typename dim>
using Config1st = push_p_config<MparticlesDouble, MfieldsC,
				InterpolateEM<Fields3d<MfieldsC::fields_t>, opt_ip_1st, dim>,
				dim, opt_ip_1st, opt_order_1st,
				CurrentNone>;

template<typename dim>
using Config1vbecDouble = push_p_config<MparticlesDouble, MfieldsC,
					InterpolateEM<Fields3d<MfieldsC::fields_t>, opt_ip_1st_ec, dim>,
					dim, opt_ip_1st_ec, opt_order_1st,
					Current1vbVar1>;

template<typename dim>
using Config1vbecSingle = push_p_config<MparticlesSingle, MfieldsSingle,
					InterpolateEM<Fields3d<MfieldsSingle::fields_t>, opt_ip_1st_ec, dim>,
					dim, opt_ip_1st_ec, opt_order_1st,
					Current1vbVar1>;

using Config1vbecSingleXZ = push_p_config<MparticlesSingle, MfieldsSingle,
					  InterpolateEM<Fields3d<MfieldsSingle::fields_t, dim_xz>, opt_ip_1st_ec, dim_xyz>,
					  dim_xyz, opt_ip_1st_ec, opt_order_1st,
					  Current1vbSplit, CacheFieldsNone,
					  dim_xz>;
using Config1vbecSingle1 = push_p_config<MparticlesSingle, MfieldsSingle,
					 InterpolateEM<Fields3d<MfieldsSingle::fields_t, dim_1>, opt_ip_1st_ec, dim_1>,
					 dim_1, opt_ip_1st_ec, opt_order_1st,
					 Current1vbVar1, CacheFieldsNone,
					 dim_1>;
