
#pragma once

// ======================================================================
// VpicHydroArrayOps

template<typename HydroArray>
struct VpicHydroArrayOps
{
  static void clear(HydroArray& ha)       { clear_hydro_array(&ha);       }
  static void synchronize(HydroArray& ha) { synchronize_hydro_array(&ha); }
};

