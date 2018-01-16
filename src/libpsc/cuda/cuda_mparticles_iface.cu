
#include "cuda_iface.h"
#include "cuda_mparticles.h"

#if 1
#define dprintf(args...) mprintf(args)
#else
#define dprintf(args...) do {} while (0)
#endif

psc_mparticles_cuda::psc_mparticles_cuda(Grid_t& grid, const Int3& bs)
{
  mprintf("CMPRTS: ctor\n");
  cmprts_ = new cuda_mparticles(grid, bs);
}

psc_mparticles_cuda::~psc_mparticles_cuda()
{
  mprintf("CMPRTS: dtor\n");
  delete cmprts_;
}

uint psc_mparticles_cuda::n_patches()
{
  mprintf("CMPRTS: n_patches\n");
  return cmprts_->n_patches;
}

void psc_mparticles_cuda::reserve_all(const uint *n_prts_by_patch)
{
  mprintf("CMPRTS: reserve_all\n");
  cmprts_->reserve_all(n_prts_by_patch);
}

void psc_mparticles_cuda::get_size_all(uint *n_prts_by_patch)
{
  mprintf("CMPRTS: get_size_all\n");
  cmprts_->get_size_all(n_prts_by_patch);
}

void psc_mparticles_cuda::resize_all(const uint *n_prts_by_patch)
{
  mprintf("CMPRTS: resize_all\n");
  cmprts_->resize_all(n_prts_by_patch);
}

uint psc_mparticles_cuda::get_n_prts()
{
  mprintf("CMPRTS: get_n_prts\n");
  return cmprts_->get_n_prts();
}

void psc_mparticles_cuda::set_particles(uint n_prts, uint off,
					void (*get_particle)(cuda_mparticles_prt *prt, int n, void *ctx),
					void *ctx)
{
  mprintf("CMPRTS: set_particles\n");
  cmprts_->set_particles(n_prts, off, get_particle, ctx);
}

void psc_mparticles_cuda::get_particles(uint n_prts, uint off,
					void (*put_particle)(cuda_mparticles_prt *, int, void *),
					void *ctx)
{
  mprintf("CMPRTS: get_particles\n");
  cmprts_->get_particles(n_prts, off, put_particle, ctx);
}

void psc_mparticles_cuda::to_device(float_4 *xi4, float_4 *pxi4,
				    uint n_prts, uint off)
{
  mprintf("CMPRTS: to_device\n");
  cmprts_->to_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::from_device(float_4 *xi4, float_4 *pxi4,
				      uint n_prts, uint off)
{
  mprintf("CMPRTS: from_device\n");
  cmprts_->from_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::setup_internals()
{
  mprintf("CMPRTS: setup_internals\n");
  cmprts_->setup_internals();
}

void psc_mparticles_cuda::inject(cuda_mparticles_prt *buf, uint *buf_n_by_patch)
{
  mprintf("CMPRTS: inject\n");
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

