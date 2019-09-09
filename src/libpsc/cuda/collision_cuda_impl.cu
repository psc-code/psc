
#include "collision_cuda_impl.hxx"

#include "cuda_collision.cuh"
#include "mparticles_cuda.hxx"
#include "balance.hxx"

// ======================================================================
// CollisionCuda

template<typename Mparticles, typename RngState>
CollisionCuda<Mparticles, RngState>::CollisionCuda(const Grid_t& grid, int interval, double nu)
  : fwd_{new CudaCollision<typename Mparticles::CudaMparticles, RngState>{interval, nu, int(1. / grid.norm.cori + .5), grid.dt}}, // FIXME nicell hack
    balance_generation_cnt_{-1}
{}

template<typename Mparticles, typename RngState>
void CollisionCuda<Mparticles, RngState>::operator()(Mparticles& mprts)
{
  if (psc_balance_generation_cnt > this->balance_generation_cnt_) {
    balance_generation_cnt_ = psc_balance_generation_cnt;
    reset(mprts);
  }
  
  (*fwd_)(*mprts.cmprts());
}

template<typename Mparticles, typename RngState>
void CollisionCuda<Mparticles, RngState>::reset(Mparticles& mprts)
{
  fwd_->reset(*mprts.cmprts());
}

template<typename Mparticles, typename RngState>
int CollisionCuda<Mparticles, RngState>::interval() const
{
  return fwd_->interval();
}

template struct CollisionCuda<MparticlesCuda<BS144>>;
template struct CollisionCuda<MparticlesCuda<BS444>>;

template struct CollisionCuda<MparticlesCuda<BS144>, RngStateFake>;
