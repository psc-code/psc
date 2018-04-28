
#include "collision_cuda_impl.hxx"

// ======================================================================
// cuda_collision

template<typename cuda_mparticles>
struct cuda_collision
{
  cuda_collision(int interval, double nu)
    : interval_{interval}, nu_{nu}
  {}
  
  void operator()(cuda_mparticles* cmprts)
  {}

private:
  int interval_;
  double nu_;
};

// ======================================================================
// CollisionCuda

template<typename BS>
CollisionCuda<BS>::CollisionCuda(MPI_Comm comm, int interval, double nu)
  : coll_{comm, interval, nu},
    fwd_{new cuda_collision<cuda_mparticles<BS>>{interval, nu}}
{}

template<typename BS>
void CollisionCuda<BS>::operator()(MparticlesCuda<BS>& _mprts)
{
  (*fwd_)(_mprts.cmprts());
  auto& mprts = _mprts.template get_as<MparticlesSingle>();
  sort_(mprts);
  coll_(mprts);
  _mprts.put_as(mprts);
}

template struct CollisionCuda<BS144>;
