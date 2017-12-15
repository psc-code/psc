
#ifndef VPIC_HYDRO_ARRAY_H
#define VPIC_HYDRO_ARRAY_H

// ======================================================================
// VpicHydroArray

template<class HydroArrayBase>
struct VpicHydroArray : HydroArrayBase
{
  typedef HydroArrayBase Base;
  using typename Base::Grid;

  static VpicHydroArray* create(Grid *grid)
  {
    return static_cast<VpicHydroArray*>(Base::create(grid));
  }
  
  // use VPIC implementations
  void clear()       { clear_hydro_array(this);       }
  void synchronize() { synchronize_hydro_array(this); }
};

#endif

