
#include "bnd_particles_cuda_impl.hxx"
#include "cuda_bndp.h"

// ----------------------------------------------------------------------
// ctor

template<typename BS, typename DIM>
BndParticlesCuda<BS, DIM>::BndParticlesCuda(struct mrc_domain *domain, const Grid_t& grid)
  : Base(domain, grid),
  cbndp_(new cuda_bndp<BS, DIM>(grid))
{}

// ----------------------------------------------------------------------
// dtor

template<typename BS, typename DIM>
BndParticlesCuda<BS, DIM>::~BndParticlesCuda()
{
  delete cbndp_;
}

// ----------------------------------------------------------------------
// reset

template<typename BS, typename DIM>
void BndParticlesCuda<BS, DIM>::reset()
{
  Base::reset();
  delete(cbndp_);
  cbndp_ = new cuda_bndp<BS, DIM>(ppsc->grid());
}

// ----------------------------------------------------------------------
// operator()

template<typename BS, typename DIM>
void BndParticlesCuda<BS, DIM>::operator()(MparticlesCuda<BS>& mprts)
{
  if (psc_balance_generation_cnt > this->balance_generation_cnt_) {
    reset();
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

template<typename BS, typename DIM>
void BndParticlesCuda<BS, DIM>::exchange_particles(MparticlesBase& mprts_base)
{
  auto& mprts = mprts_base.get_as<MparticlesCuda<BS>>();
  (*this)(mprts);
  mprts_base.put_as(mprts);
}

template struct BndParticlesCuda<BS144, dim_yz>;
template struct BndParticlesCuda<BS444, dim_xyz>;