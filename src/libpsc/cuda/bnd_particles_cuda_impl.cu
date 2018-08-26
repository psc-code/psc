
#include "bnd_particles_cuda_impl.hxx"
#include "cuda_bndp.h"

// ----------------------------------------------------------------------
// ctor

template<typename Mparticles, typename DIM>
BndParticlesCuda<Mparticles, DIM>::BndParticlesCuda(const Grid_t& grid)
  : Base(grid),
    cbndp_(new cuda_bndp<typename Mparticles::CudaMparticles, DIM>(grid))
{}

// ----------------------------------------------------------------------
// dtor

template<typename Mparticles, typename DIM>
BndParticlesCuda<Mparticles, DIM>::~BndParticlesCuda()
{
  delete cbndp_;
}

// ----------------------------------------------------------------------
// reset

template<typename Mparticles, typename DIM>
void BndParticlesCuda<Mparticles, DIM>::reset(const Grid_t& grid)
{
  Base::reset(grid);
  delete(cbndp_);
  cbndp_ = new cuda_bndp<typename Mparticles::CudaMparticles, DIM>(grid);
}

// ----------------------------------------------------------------------
// operator()

template<typename Mparticles, typename DIM>
void BndParticlesCuda<Mparticles, DIM>::operator()(Mparticles& mprts)
{
  if (psc_balance_generation_cnt > this->balance_generation_cnt_) {
    reset(mprts.grid());
  }
  
  static int pr_A, pr_B;
  if (!pr_A) {
    pr_A = prof_register("xchg_mprts_prep", 1., 0, 0);
    pr_B = prof_register("xchg_mprts_post", 1., 0, 0);
  }
  
  prof_restart(pr_time_step_no_comm);
  prof_start(pr_A);
  cbndp_->prep(mprts.cmprts());
  prof_stop(pr_A);
  
  this->process_and_exchange(mprts, cbndp_->bufs_);
  
  prof_restart(pr_time_step_no_comm);
  prof_start(pr_B);
  cbndp_->post(mprts.cmprts());
  prof_stop(pr_B);
  prof_stop(pr_time_step_no_comm);
}

// ----------------------------------------------------------------------
// exchange_particles

template<typename Mparticles, typename DIM>
void BndParticlesCuda<Mparticles, DIM>::exchange_particles(MparticlesBase& mprts_base)
{
  auto& mprts = mprts_base.get_as<Mparticles>();
  (*this)(mprts);
  mprts_base.put_as(mprts);
}

template struct BndParticlesCuda<MparticlesCuda<BS144>, dim_yz>;
template struct BndParticlesCuda<MparticlesCuda<BS444>, dim_xyz>;
