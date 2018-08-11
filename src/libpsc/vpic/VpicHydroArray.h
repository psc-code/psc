
#pragma once

// ======================================================================
// VpicHydroArrayOps

template<typename HydroArray>
struct VpicHydroArrayOps
{
  static void clear(HydroArray& ha)       { clear_hydro_array(&ha);       }
  static void synchronize(HydroArray& ha) { synchronize_hydro_array(&ha); }
};
  
// ======================================================================
// VpicHydroArray

template<class HydroArrayBase>
struct VpicHydroArray : HydroArrayBase
{
  using Self = VpicHydroArray<HydroArrayBase>;
  typedef HydroArrayBase Base;
  using typename Base::Grid;

  using Base::Base;
  
  void clear() { return VpicHydroArrayOps<Self>::clear(*this); }
  void synchronize() { return VpicHydroArrayOps<Self>::synchronize(*this); }
};

