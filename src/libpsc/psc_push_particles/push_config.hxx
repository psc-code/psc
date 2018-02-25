
#include "dim.hxx"

template<typename MP, typename MF, typename DIM>
struct Config
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using dim = DIM;
};

#include "psc_particles_double.h"
#include "psc_fields_c.h"

using Config2ndXYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xyz>;
using Config2ndXY = Config<PscMparticlesDouble, PscMfieldsC, dim_xy>;
using Config2ndXZ = Config<PscMparticlesDouble, PscMfieldsC, dim_xz>;
using Config2ndYZ = Config<PscMparticlesDouble, PscMfieldsC, dim_yz>;


