
#include "dim.hxx"
#include "inc_defs.h"

struct CacheFieldsNone;
struct CacheFields;

template<typename MP, typename MF, typename DIM,
	 typename ORDER = opt_order_2nd,
	 typename CF = CacheFieldsNone>
struct Config
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using dim = DIM;
  using order = ORDER;
  using CacheFields = CF;
};

#include "psc_particles_double.h"
#include "psc_fields_c.h"

using Config2ndXYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xyz>;
using Config2ndXY = Config<PscMparticlesDouble, PscMfieldsC, dim_xy>;
using Config2ndXZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xz>;
using Config2ndYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz>;
using Config2ndY = Config<PscMparticlesDouble, PscMfieldsC, dim_y>;
using Config2ndZ = Config<PscMparticlesDouble, PscMfieldsC, dim_z>;

using Config2ndDoubleYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_order_2nd, CacheFields>;

using Config1stXZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xz, opt_order_1st>;
using Config1stYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz, opt_order_1st>;



