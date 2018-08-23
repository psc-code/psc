
#ifndef PSC_SIMULATION_BASE_H
#define PSC_SIMULATION_BASE_H

template<class Particles, class MaterialList>
class PscSimulationMixin
{
  typedef typename Particles::Grid Grid;

public:
  PscSimulationMixin()
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

};


#endif
