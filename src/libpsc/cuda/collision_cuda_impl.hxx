
#pragma once

#include "../libpsc/psc_sort/psc_sort_impl.hxx"

template<typename BS>
class CollisionCuda : CollisionBase
{
public:
  CollisionCuda(MPI_Comm comm, int interval, double nu)
    : coll_{comm, interval, nu}
  {}
  
  virtual void run(PscMparticlesBase mprts_base) { assert(0); }

  void operator()(MparticlesCuda<BS>& _mprts)
  {
    auto& mprts = _mprts.template get_as<MparticlesSingle>();
    sort_(mprts);
    coll_(mprts);
    _mprts.put_as(mprts);
  }

private:
  SortCountsort2<MparticlesSingle> sort_;
  Collision_<MparticlesSingle, MfieldsSingle> coll_;
};

