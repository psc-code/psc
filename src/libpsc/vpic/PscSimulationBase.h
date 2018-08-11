
#ifndef PSC_SIMULATION_BASE_H
#define PSC_SIMULATION_BASE_H

template<class Particles, class MaterialList>
class PscSimulationMixin
{
  typedef typename Particles::Grid Grid;

public:
  PscSimulationMixin()
  {
    grid_ = Grid::create();
  }
  
  Grid*& getGrid()                    { return grid_; }
  
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
};


#endif
