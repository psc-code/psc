
#ifndef PSC_SIMULATION_BASE_H
#define PSC_SIMULATION_BASE_H

#include "VpicDiag.h"

template<class Particles>
class PscSimulationMixin : protected vpic_simulation
{
  typedef typename Particles::Grid Grid;
  typedef typename Particles::FieldArray FieldArray;
  typedef typename Particles::Interpolator Interpolator;
  typedef typename Particles::Accumulator Accumulator;
  typedef typename Particles::HydroArray HydroArray;
  typedef typename FieldArray::MaterialList MaterialList;

public:
  Grid*& getGrid()                 { return *reinterpret_cast<Grid **>(&grid); }
  MaterialList& getMaterialList()  { return material_list_; }
  FieldArray*& getFieldArray()     { return field_array_; }
  Interpolator*& getInterpolator() { return interpolator_; }
  Accumulator*& getAccumulator()   { return accumulator_; }
  HydroArray*& getHydroArray()     { return hydro_array_;  }
  Particles& getParticles()        { return particles_; }
  
  void getParams(int& num_step_,
		 int& clean_div_e_interval_,
		 int& clean_div_b_interval_,
		 int& sync_shared_interval_,
		 int& num_div_e_round_,
		 int& num_div_b_round_,
		 int& status_interval_)
  {
    assert(0);
  }

  void setParams(int num_step_, int status_interval_,
		 int sync_shared_interval_, int clean_div_e_interval_,
		 int clean_div_b_interval_)
  {
  }

  void setTopology(int px_, int py_, int pz_)
  {
    px = px_; py = py_; pz = pz_;
  }

  //FIXME, those should probably be in a separate mixin...
  void initialization(int argc, char **argv)
  {
  }

  void diagnostics()
  {
  }
  
  void emitter()
  {
  }

  void collision_run()
  {
  }

  void current_injection()
  {
  }

  void field_injection()
  {
  }

private:
  MaterialList material_list_;
  FieldArray *field_array_;
  Interpolator *interpolator_;
  Accumulator *accumulator_;
  HydroArray *hydro_array_;
  Particles particles_;
};

template<class Particles>
class PscSimulationBase
{
};


#endif
