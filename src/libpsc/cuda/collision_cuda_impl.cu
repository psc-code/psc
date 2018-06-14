
#include "collision_cuda_impl.hxx"

#include "cuda_collision.cuh"

// ======================================================================
// CollisionCuda

template<typename BS, typename RngState>
CollisionCuda<BS, RngState>::CollisionCuda(MPI_Comm comm, int interval, double nu)
  : fwd_{new CudaCollision<cuda_mparticles<BS>, RngState>{interval, nu, ppsc->prm.nicell, ppsc->dt}}
{}

template<typename BS, typename RngState>
void CollisionCuda<BS, RngState>::operator()(MparticlesCuda<BS>& _mprts)
{
  (*fwd_)(*_mprts.cmprts());
}

template struct CollisionCuda<BS144>;
template struct CollisionCuda<BS444>;

template struct CollisionCuda<BS144, RngStateFake>;
