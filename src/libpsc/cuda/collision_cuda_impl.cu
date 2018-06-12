
#include "collision_cuda_impl.hxx"

#include "cuda_collision.cuh"

// ======================================================================
// CollisionCuda

template<typename BS>
CollisionCuda<BS>::CollisionCuda(MPI_Comm comm, int interval, double nu)
  : fwd_{new cuda_collision<cuda_mparticles<BS>>{interval, nu, ppsc->prm.nicell, ppsc->dt}}
{}

template<typename BS>
void CollisionCuda<BS>::operator()(MparticlesCuda<BS>& _mprts)
{
  (*fwd_)(*_mprts.cmprts());
}

template struct CollisionCuda<BS144>;
template struct CollisionCuda<BS444>;
