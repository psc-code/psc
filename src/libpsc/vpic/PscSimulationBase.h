
#ifndef PSC_SIMULATION_BASE_H
#define PSC_SIMULATION_BASE_H

#include "VpicDiag.h"

template<class FieldArray, class Particles, class Interpolator, class Accumulator,
	 class HydroArray>
class PscSimulationBase : protected vpic_simulation
{
public:
  PscSimulationBase()
  {
  }

  Grid*& getGrid()
  {
    return *reinterpret_cast<Grid **>(&grid);
  }

  MaterialList& getMaterialList()
  {
    return material_list_;
  }

  FieldArray*& getFieldArray()
  {
    return field_array_;
  }

  Interpolator*& getInterpolator()
  {
    return interpolator_;
  }
  
  Accumulator*& getAccumulator()
  {
    return accumulator_;
  }

  HydroArray*& getHydroArray()
  {
    return hydro_array_;
  }

  VpicParticles& getParticles()
  {
    return *reinterpret_cast<VpicParticles *>(&species_list);
  }
  
  void emitter()
  {
  }

  void collision_run()
  {
  }

  void user_current_injection()
  {
  }

  void user_field_injection()
  {
  }
  
  void setParams(int num_step_, int status_interval_,
		 int sync_shared_interval_, int clean_div_e_interval_,
		 int clean_div_b_interval_)
  {
  }

  void setTopology(int px_, int py_, int pz_)
  {
  }

 private:
  MaterialList material_list_;
  FieldArray *field_array_;
  Interpolator *interpolator_;
  Accumulator *accumulator_;
  HydroArray *hydro_array_;
};


#endif
