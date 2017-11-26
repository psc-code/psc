
#ifndef VPIC_HYDRO_ARRAY_H
#define VPIC_HYDRO_ARRAY_H

// ======================================================================
// VpicHydroArray

template<class HydroArrayBase>
struct VpicHydroArray : HydroArrayBase
{
  typedef HydroArrayBase Base;

  using Base::Base;

  // use VPIC implementations
  void clear()       { clear_hydro_array(this);       }
  void synchronize() { synchronize_hydro_array(this); }
};

#endif

