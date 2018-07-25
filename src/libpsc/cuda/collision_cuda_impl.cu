
#include "collision_cuda_impl.hxx"

#include "cuda_collision.cuh"

// ======================================================================
// CollisionCuda

template<typename Mparticles, typename RngState>
CollisionCuda<Mparticles, RngState>::CollisionCuda(MPI_Comm comm, int interval, double nu)
  : fwd_{new CudaCollision<typename Mparticles::CudaMparticles, RngState>{interval, nu, int(1. / ppsc->grid().norm.cori + .5), ppsc->grid().dt}} // FIXME nicell hack
{}

template<typename Mparticles, typename RngState>
void CollisionCuda<Mparticles, RngState>::operator()(Mparticles& _mprts)
{
  (*fwd_)(*_mprts.cmprts());
}

template struct CollisionCuda<MparticlesCuda<BS144>>;
template struct CollisionCuda<MparticlesCuda<BS444>>;

template struct CollisionCuda<MparticlesCuda<BS144>, RngStateFake>;
