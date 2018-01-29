
#ifndef CUDA_BNDP_IFACE_H
#define CUDA_BNDP_IFACE_H

#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_bnd_particles_private.h"

#include "bnd_particles_impl.hxx"
#include "cuda_bndp.h"

struct psc_bnd_particles_cuda;

template<typename MP>
struct bnd_particles_policy_cuda
{
  using mparticles_t = MP;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  // ----------------------------------------------------------------------
  // ctor
  
  bnd_particles_policy_cuda()
    : cbndp_(new cuda_bndp())
  {}

  // ----------------------------------------------------------------------
  // setup

  void setup(Grid_t& grid)
  {
    cbndp_->setup(grid);
  }

protected:
  // ----------------------------------------------------------------------
  // exchange_mprts_prep
  
  void exchange_mprts_prep(ddcp_t* ddcp, mparticles_t mprts)
  {
    cbndp_->prep(ddcp, mprts->cmprts_);
  }

  // ----------------------------------------------------------------------
  // exchange_mprts_post
  
  void exchange_mprts_post(ddcp_t* ddcp, mparticles_t mprts)
  {
    cbndp_->post(ddcp, mprts->cmprts_);
  }

private:
  std::unique_ptr<cuda_bndp> cbndp_;
};

struct psc_bnd_particles_cuda : psc_bnd_particles_sub<mparticles_cuda_t,
						      bnd_particles_policy_cuda<mparticles_cuda_t>>
{
  using Base = psc_bnd_particles_sub<mparticles_cuda_t,
				     bnd_particles_policy_cuda<mparticles_cuda_t>>;

  using Base::setup;
  
  // ----------------------------------------------------------------------
  // exchange_particles

  void exchange_particles(mparticles_cuda_t mprts)
  {
    Base::exchange_particles(mprts);
  }

  // ======================================================================
  // interface to psc_bnd_particles
  // repeated here since there's no way to do this somehow virtual at
  // this spoint
  
  // ----------------------------------------------------------------------
  // create

  static void create(struct psc_bnd_particles *bnd)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);

    new(sub) psc_bnd_particles_cuda();
  }
  
  // ----------------------------------------------------------------------
  // setup
  
  static void setup(struct psc_bnd_particles *bnd)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);
    
    sub->setup(bnd->psc->mrc_domain);
    static_cast<policy_t*>(sub)->setup(bnd->psc->grid);
  }

  // ----------------------------------------------------------------------
  // destroy

  static void destroy(struct psc_bnd_particles *bnd)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);

    sub->~psc_bnd_particles_cuda();
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
