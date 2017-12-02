
#ifndef PSC_SIMULATION_BASE_H
#define PSC_SIMULATION_BASE_H

#include "VpicDiag.h"

template<class Particles>
class PscSimulationMixin
{
  typedef typename Particles::Grid Grid;
  typedef typename Particles::FieldArray FieldArray;
  typedef typename Particles::Interpolator Interpolator;
  typedef typename Particles::Accumulator Accumulator;
  typedef typename Particles::HydroArray HydroArray;
  typedef typename Particles::ParticleBcList ParticleBcList;
  typedef typename FieldArray::MaterialList MaterialList;

public:
  PscSimulationMixin()
  {
    grid_ = Grid::create();
  }
  
  Grid*& getGrid()                    { return grid_; }
  MaterialList& getMaterialList()     { return material_list_; }
  FieldArray*& getFieldArray()        { return field_array_; }
  Interpolator*& getInterpolator()    { return interpolator_; }
  Accumulator*& getAccumulator()      { return accumulator_; }
  HydroArray*& getHydroArray()        { return hydro_array_;  }
  Particles& getParticles()           { return particles_; }
  ParticleBcList& getParticleBcList() { return particle_bc_list_; }
  
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
  Grid* grid_;
  MaterialList material_list_;
  FieldArray *field_array_;
  Interpolator *interpolator_;
  Accumulator *accumulator_;
  HydroArray *hydro_array_;
  Particles particles_;
  ParticleBcList particle_bc_list_;
};

template<class Particles>
class PscSimulationBase
{
};


#endif
