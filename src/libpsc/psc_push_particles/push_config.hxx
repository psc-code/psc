
#pragma once

#include "dim.hxx"
#include "inc_defs.h"

template<typename dim>
struct CacheFieldsNone;

template<typename dim>
struct CacheFields;

template<typename MP, typename MF, typename D,
	 typename IP, typename O,
	 typename CALCJ = opt_calcj_esirkepov,
	 typename OPT_EXT = opt_ext_none,
	 typename CF = CacheFieldsNone<D>>
struct push_p_config
{
  using Mparticles = typename MP::sub_t;
  using Mfields = typename MF::sub_t;
  using dim = D;
  using ip = IP;
  using order = O;
  using calcj = CALCJ;
  using ext = OPT_EXT;
  using CacheFields = CF;
};

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

template<typename dim>
using Config2nd = push_p_config<PscMparticlesDouble, PscMfieldsC, dim, opt_ip_2nd, opt_order_2nd>;

using Config2ndDoubleYZ = push_p_config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_ip_2nd, opt_order_2nd,
					opt_calcj_esirkepov, opt_ext_none, CacheFields<dim_yz>>;

template<typename dim>
using Config1st = push_p_config<PscMparticlesDouble, PscMfieldsC, dim, opt_ip_1st, opt_order_1st>;

template<typename dim>
using Config1vbecDouble = push_p_config<PscMparticlesDouble, PscMfieldsC, dim, opt_ip_1st_ec,
					opt_order_1st, opt_calcj_1vb_var1>;

template<typename dim>
using Config1vbecSingle = push_p_config<PscMparticlesSingle, PscMfieldsSingle, dim, opt_ip_1st_ec,
					opt_order_1st, opt_calcj_1vb_var1>;

