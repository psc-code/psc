
#include "cuda_iface.h"
#include "cuda_mparticles.h"

psc_mparticles_cuda::psc_mparticles_cuda(Grid_t& grid, const Int3& bs)
{
  cmprts_ = new cuda_mparticles(grid, bs);
}

psc_mparticles_cuda::~psc_mparticles_cuda()
{
  delete cmprts_;
}

uint psc_mparticles_cuda::n_patches()
{
  return cmprts_->n_patches;
}

void psc_mparticles_cuda::reserve_all(const uint *n_prts_by_patch)
{
  cmprts_->reserve_all(n_prts_by_patch);
}

void psc_mparticles_cuda::get_size_all(uint *n_prts_by_patch)
{
  cmprts_->get_size_all(n_prts_by_patch);
}

void psc_mparticles_cuda::resize_all(const uint *n_prts_by_patch)
{
  cmprts_->resize_all(n_prts_by_patch);
}

uint psc_mparticles_cuda::get_n_prts()
{
  return cmprts_->get_n_prts();
}

void psc_mparticles_cuda::set_particles(uint n_prts, uint off,
					void (*get_particle)(cuda_mparticles_prt *prt, int n, void *ctx),
					void *ctx)
{
  cmprts_->set_particles(n_prts, off, get_particle, ctx);
}

void psc_mparticles_cuda::get_particles(uint n_prts, uint off,
					void (*put_particle)(cuda_mparticles_prt *, int, void *),
					void *ctx)
{
  cmprts_->get_particles(n_prts, off, put_particle, ctx);
}

void psc_mparticles_cuda::to_device(float_4 *xi4, float_4 *pxi4,
				    uint n_prts, uint off)
{
  cmprts_->to_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::from_device(float_4 *xi4, float_4 *pxi4,
				      uint n_prts, uint off)
{
  cmprts_->from_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::setup_internals()
{
  cmprts_->setup_internals();
}

void psc_mparticles_cuda::inject(cuda_mparticles_prt *buf, uint *buf_n_by_patch)
{
  cmprts_->inject(buf, buf_n_by_patch);
}

const particle_cuda_real_t* psc_mparticles_cuda::patch_get_b_dxi(int p)
{
  return cmprts_->patch_get_b_dxi(p);
}

const int* psc_mparticles_cuda::patch_get_b_mx(int p)
{
  return cmprts_->patch_get_b_mx(p);
}

psc_particle_cuda_buf_t* psc_mparticles_cuda::bnd_get_buffer(int p)
{
  return cmprts_->bnd_get_buffer(p);
}

void psc_mparticles_cuda::bnd_prep()
{
  cmprts_->bnd_prep();
}

void psc_mparticles_cuda::bnd_post()
{
  cmprts_->bnd_post();
}

