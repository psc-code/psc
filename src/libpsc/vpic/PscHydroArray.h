
#ifndef PSC_HYDRO_ARRAY_H
#define PSC_HYDRO_ARRAY_H

// ======================================================================
// PscHydroArray

template<class HydroArrayBase>
struct PscHydroArray : HydroArrayBase
{
  typedef HydroArrayBase Base;

  using Base::Base;

  void clear()       { clear_hydro_array(this);       }
  void synchronize() { synchronize_hydro_array(this); }
};

#endif

