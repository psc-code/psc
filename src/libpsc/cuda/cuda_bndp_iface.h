
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

protected:
  // ----------------------------------------------------------------------
  // exchange_mprts_prep
  
  void exchange_mprts_prep(ddcp_t* ddcp, mparticles_t mprts)
  {
    if (!is_setup_) {
      cbndp_->setup(ppsc->grid);
      is_setup_ = true;
    }
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
  bool is_setup_ = false; // FIXME, setup should be done at ctor time
};

struct psc_bnd_particles_cuda : psc_bnd_particles_sub<mparticles_cuda_t,
						      bnd_particles_policy_cuda<mparticles_cuda_t>>
{
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
  // destroy

  static void destroy(struct psc_bnd_particles *bnd)
  {
    auto sub = static_cast<psc_bnd_particles_cuda*>(bnd->obj.subctx);

    sub->~psc_bnd_particles_cuda();
  }
  
};

#endif
