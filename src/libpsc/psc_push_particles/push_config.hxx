
#pragma once

#include "dim.hxx"
#include "inc_defs.h"

struct CacheFieldsNone;
struct CacheFields;

template<typename MP, typename MF, typename D,
	 typename IP, typename O,
	 typename CALCJ, typename OPT_EXT = opt_ext_none>
struct push_p_config
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using dim = D;
  using ip = IP;
  using order = O;
  using calcj = CALCJ;
  using ext = OPT_EXT;
};

template<typename MP, typename MF, typename D,
	 typename IP, typename O = opt_order_2nd,
	 typename CF = CacheFieldsNone>
struct Config
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using dim = D;
  using ip = IP;
  using order = O;
  using CacheFields = CF;
};

#include "psc_particles_double.h"
#include "psc_fields_c.h"

using Config2ndXYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xyz, opt_ip_2nd>;
using Config2ndXY = Config<PscMparticlesDouble, PscMfieldsC, dim_xy, opt_ip_2nd>;
using Config2ndXZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xz, opt_ip_2nd>;
using Config2ndYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_ip_2nd>;
using Config2ndY = Config<PscMparticlesDouble, PscMfieldsC, dim_y, opt_ip_2nd>;
using Config2ndZ = Config<PscMparticlesDouble, PscMfieldsC, dim_z, opt_ip_2nd>;
using Config2nd1 = Config<PscMparticlesDouble, PscMfieldsC, dim_1, opt_ip_2nd>;

using Config2ndDoubleYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_ip_2nd, opt_order_2nd, CacheFields>;

using Config1stXZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xz, opt_ip_1st, opt_order_1st>;
using Config1stYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_ip_1st, opt_order_1st>;
