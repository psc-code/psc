
#pragma once

#include "vpic/vpic.h"
#include "psc_vpic_bits.h"

template<typename _Mparticles, typename _MfieldsHydro, typename _MfieldsInterpolator>
struct VpicHydroOps
{
  using Mparticles = _Mparticles;
  using MfieldsHydro = _MfieldsHydro;
  using MfieldsInterpolator = _MfieldsInterpolator;
  
  static void clear(MfieldsHydro& hydro)       { ::clear_hydro_array(hydro); }
  static void synchronize(MfieldsHydro& hydro) { ::synchronize_hydro_array(hydro); }
  static void accumulate_hydro_p(MfieldsHydro& hydro, typename Mparticles::const_iterator sp,
				 /*const*/ MfieldsInterpolator& interpolator)
  {
    ::accumulate_hydro_p(hydro, &*sp, interpolator.getPatch(0).ip());
  }
};
