
#pragma once

// ======================================================================
// VpicHydroArrayOps

template<typename HydroArray, typename MfieldsHydro>
struct VpicHydroArrayOps
{
  static void clear(MfieldsHydro& hydro)       { clear_hydro_array(&hydro.vhydro()); }
  static void synchronize(MfieldsHydro& hydro) { synchronize_hydro_array(&hydro.vhydro()); }
};

