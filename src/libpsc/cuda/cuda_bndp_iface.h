
#ifndef CUDA_BNDP_IFACE_H
#define CUDA_BNDP_IFACE_H

#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_bnd_particles_private.h"

#include "bnd_particles_impl.hxx"
#include "cuda_bndp.h"

extern int pr_time_step_no_comm;

struct psc_bnd_particles_cuda : psc_bnd_particles_sub<PscMparticlesCuda>
{
  using Base = psc_bnd_particles_sub<PscMparticlesCuda>;

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

  void exchange_particles(PscMparticlesCuda mprts)
  {
    static int pr_A, pr_B;
    if (!pr_A) {
      pr_A = prof_register("xchg_mprts_prep", 1., 0, 0);
      pr_B = prof_register("xchg_mprts_post", 1., 0, 0);
    }
    
    prof_restart(pr_time_step_no_comm);
    prof_start(pr_A);
    cbndp_->prep(ddcp, mprts->cmprts());
    prof_stop(pr_A);
    
    process_and_exchange(mprts);
    
    prof_restart(pr_time_step_no_comm);
    prof_start(pr_B);
    cbndp_->post(ddcp, mprts->cmprts());
    prof_stop(pr_B);
    prof_stop(pr_time_step_no_comm);
  }

private:
  std::unique_ptr<cuda_bndp> cbndp_;

public:
  // ======================================================================
  // interface to psc_bnd_particles
  // repeated here since there's no way to do this somehow virtual at
  // this spoint
  
  // ----------------------------------------------------------------------
  // reset
  
  static void reset(struct psc_bnd_particles *bnd)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);
    
    sub->reset(bnd->psc->mrc_domain, bnd->psc->grid());
  }

  // ----------------------------------------------------------------------
  // exchange_particles

  static void exchange_particles(struct psc_bnd_particles *bnd,
				 struct psc_mparticles *mprts_base)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);
    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
    
    sub->exchange_particles(mprts);
    
    mprts.put_as(mprts_base);
  }
};

#endif
