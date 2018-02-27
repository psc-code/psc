
#pragma once

#include "dim.hxx"
#include "inc_defs.h"

template<typename fields_t, typename dim>
struct CacheFieldsNone;

template<typename fields_t, typename dim>
struct CacheFields;

#include "fields.hxx"

template<typename MP, typename MF, typename D,
	 typename IP, typename O,
	 typename CALCJ = opt_calcj_esirkepov,
	 typename OPT_EXT = opt_ext_none,
	 typename CF = CacheFieldsNone<typename MF::fields_t, D>,
	 typename EM = Fields3d<typename MF::fields_t>>
struct push_p_config
{
  using Mparticles = MP;
  using Mfields = MF;
  using dim = D;
  using ip = IP;
  using order = O;
  using calcj = CALCJ;
  using ext = OPT_EXT;
  using CacheFields = CF;

  using FieldsEM = EM;
};

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

template<typename dim>
using Config2nd = push_p_config<MparticlesDouble, MfieldsC, dim, opt_ip_2nd, opt_order_2nd>;

using Config2ndDoubleYZ = push_p_config<MparticlesDouble, MfieldsC, dim_yz, opt_ip_2nd, opt_order_2nd,
					opt_calcj_esirkepov, opt_ext_none, CacheFields<MfieldsC::fields_t, dim_yz>>;

template<typename dim>
using Config1st = push_p_config<MparticlesDouble, MfieldsC, dim, opt_ip_1st, opt_order_1st>;

template<typename dim>
using Config1vbecDouble = push_p_config<MparticlesDouble, MfieldsC, dim, opt_ip_1st_ec,
					opt_order_1st, opt_calcj_1vb_var1>;

template<typename dim>
using Config1vbecSingle = push_p_config<MparticlesSingle, MfieldsSingle, dim, opt_ip_1st_ec,
					opt_order_1st, opt_calcj_1vb_var1>;

using Config1vbecSingleXZ = push_p_config<MparticlesSingle, MfieldsSingle, dim_xyz, opt_ip_1st_ec, opt_order_1st,
					  opt_calcj_1vb_split, opt_ext_none, CacheFieldsNone<MfieldsC::fields_t, dim_xyz>,
					  Fields3d<typename MfieldsSingle::fields_t, dim_xz>>;
using Config1vbecSingle1 = push_p_config<MparticlesSingle, MfieldsSingle, dim_1, opt_ip_1st_ec, opt_order_1st,
					 opt_calcj_1vb_var1, opt_ext_none, CacheFieldsNone<MfieldsC::fields_t, dim_1>,
					 Fields3d<typename MfieldsSingle::fields_t, dim_1>>;
