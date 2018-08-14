
#pragma once

// ======================================================================
// VpicHydroArrayOps

template<typename MfieldsHydro>
struct VpicHydroArrayOps
{
  static void clear(MfieldsHydro& hydro)       { clear_hydro_array(hydro); }
  static void synchronize(MfieldsHydro& hydro) { synchronize_hydro_array(hydro); }
};

