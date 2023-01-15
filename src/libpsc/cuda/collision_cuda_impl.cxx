
#include "collision_cuda_impl.hxx"

#include "balance.hxx"
#include "cuda_collision.hxx"
#include "mparticles_cuda.hxx"

// ======================================================================
// CollisionCuda

template <typename Mparticles, typename RngState>
CollisionCuda<Mparticles, RngState>::CollisionCuda(const Grid_t& grid,
                                                   int interval, double nu)
  : fwd_{new CudaCollision<typename Mparticles::CudaMparticles, RngState>{
      interval, nu, int(1. / grid.norm.cori + .5), grid.dt}}
// FIXME nicell hack
{}

template <typename Mparticles, typename RngState>
void CollisionCuda<Mparticles, RngState>::operator()(Mparticles& mprts)
{
  (*fwd_)(*mprts.cmprts());
}

template <typename Mparticles, typename RngState>
int CollisionCuda<Mparticles, RngState>::interval() const
{
  return fwd_->interval();
}

template struct CollisionCuda<MparticlesCuda<BS144>>;
template struct CollisionCuda<MparticlesCuda<BS444>>;

// template struct CollisionCuda<MparticlesCuda<BS144>, RngStateFake>;
