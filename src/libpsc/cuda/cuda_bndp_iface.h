
#ifndef CUDA_BNDP_IFACE_H
#define CUDA_BNDP_IFACE_H

#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_bnd_particles_private.h"

#include "bnd_particles_impl.hxx"
#include "cuda_bndp.h"

extern int pr_time_step_no_comm;

struct psc_bnd_particles_cuda : psc_bnd_particles_sub<MparticlesCuda>
{
  using Base = psc_bnd_particles_sub<MparticlesCuda>;
  using mparticles_t = PscMparticlesCuda;

  // ----------------------------------------------------------------------
  // ctor
  
  psc_bnd_particles_cuda(struct mrc_domain *domain, const Grid_t& grid)
    : Base(domain, grid),
      cbndp_(new cuda_bndp(grid))
  {}

  // ----------------------------------------------------------------------
  // reset
  
  void reset(struct mrc_domain *domain, const Grid_t& grid) override
  {
    Base::reset(domain, grid);
    //cbndp_->setup(grid);
  }

  // ----------------------------------------------------------------------
  // exchange_particles

  void exchange_particles(PscMparticlesBase mprts_base) override
  {
    static int pr_A, pr_B;
    if (!pr_A) {
      pr_A = prof_register("xchg_mprts_prep", 1., 0, 0);
      pr_B = prof_register("xchg_mprts_post", 1., 0, 0);
    }
    
    auto mprts = mprts_base.get_as<mparticles_t>();

    prof_restart(pr_time_step_no_comm);
    prof_start(pr_A);
    cbndp_->prep(ddcp, mprts->cmprts());
    prof_stop(pr_A);
    
    process_and_exchange(*mprts.sub());
    
    prof_restart(pr_time_step_no_comm);
    prof_start(pr_B);
    cbndp_->post(ddcp, mprts->cmprts());
    prof_stop(pr_B);
    prof_stop(pr_time_step_no_comm);

    mprts.put_as(mprts_base);
}

private:
  std::unique_ptr<cuda_bndp> cbndp_;
};

#endif
