
#include "bnd_particles_cuda_impl.hxx"
#include "cuda_bndp.h"

// ----------------------------------------------------------------------
// ctor

BndParticlesCuda::BndParticlesCuda(struct mrc_domain *domain, const Grid_t& grid)
  : Base(domain, grid),
  cbndp_(new cuda_bndp(grid))
{}

// ----------------------------------------------------------------------
// dtor

BndParticlesCuda::~BndParticlesCuda()
{
  delete cbndp_;
}

// ----------------------------------------------------------------------
// reset

void BndParticlesCuda::reset()
{
  Base::reset();
  //cbndp_->setup(grid);
}

// ----------------------------------------------------------------------
// exchange_particles

void BndParticlesCuda::exchange_particles(PscMparticlesBase mprts_base)
{
  // FIXME, needs to check for rebalanced domain
  
  static int pr_A, pr_B;
  if (!pr_A) {
    pr_A = prof_register("xchg_mprts_prep", 1., 0, 0);
    pr_B = prof_register("xchg_mprts_post", 1., 0, 0);
  }
  
  auto& mprts = mprts_base->get_as<MparticlesCuda>();
  
  prof_restart(pr_time_step_no_comm);
  prof_start(pr_A);
  cbndp_->prep(ddcp, mprts.cmprts());
  prof_stop(pr_A);
  
  process_and_exchange(mprts);
  
  prof_restart(pr_time_step_no_comm);
  prof_start(pr_B);
  cbndp_->post(ddcp, mprts.cmprts());
  prof_stop(pr_B);
  prof_stop(pr_time_step_no_comm);
  
  mprts_base->put_as(mprts);
}

