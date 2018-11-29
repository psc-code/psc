
#pragma once

// ======================================================================
// VpicHydroArrayOps

template<typename _MfieldsHydro>
struct VpicHydroArrayOps
{
  using MfieldsHydro = _MfieldsHydro;
  
  static void clear(MfieldsHydro& hydro)       { clear_hydro_array(hydro); }
  static void synchronize(MfieldsHydro& hydro) { synchronize_hydro_array(hydro); }
};

