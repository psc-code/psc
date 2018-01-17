
#include "cuda_iface.h"
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"

#if 1
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

psc_mparticles_cuda::psc_mparticles_cuda(Grid_t& grid, const Int3& bs)
{
  dprintf("CMPRTS: ctor\n");
  cmprts_ = new cuda_mparticles(grid, bs);
}

psc_mparticles_cuda::~psc_mparticles_cuda()
{
  dprintf("CMPRTS: dtor\n");
  delete cmprts_;
}

uint psc_mparticles_cuda::n_patches()
{
  dprintf("CMPRTS: n_patches\n");
  return cmprts_->n_patches;
}

void psc_mparticles_cuda::reserve_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: reserve_all\n");
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  cmprts_->reserve_all(n_prts_by_patch);
}

void psc_mparticles_cuda::get_size_all(uint *n_prts_by_patch)
{
  dprintf("CMPRTS: get_size_all\n");
  cmprts_->get_size_all(n_prts_by_patch);
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
}

void psc_mparticles_cuda::resize_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: resize_all\n");
  cmprts_->resize_all(n_prts_by_patch);
}

uint psc_mparticles_cuda::get_n_prts()
{
  dprintf("CMPRTS: get_n_prts\n");
  return cmprts_->get_n_prts();
}

void psc_mparticles_cuda::set_particles(uint n_prts, uint off,
					void (*get_particle)(cuda_mparticles_prt *prt, int n, void *ctx),
					void *ctx)
{
  dprintf("CMPRTS: set_particles\n");
  cmprts_->set_particles(n_prts, off, get_particle, ctx);
}

void psc_mparticles_cuda::get_particles(uint n_prts, uint off,
					void (*put_particle)(cuda_mparticles_prt *, int, void *),
					void *ctx)
{
  dprintf("CMPRTS: get_particles\n");
  cmprts_->get_particles(n_prts, off, put_particle, ctx);
}

void psc_mparticles_cuda::to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off)
{
  dprintf("CMPRTS: to_device\n");
  cmprts_->to_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off)
{
  dprintf("CMPRTS: from_device\n");
  cmprts_->from_device(xi4, pxi4, n_prts, off);
}

void psc_mparticles_cuda::setup_internals()
{
  dprintf("CMPRTS: setup_internals\n");
  cmprts_->setup_internals();
}

void psc_mparticles_cuda::inject(cuda_mparticles_prt *buf, uint *buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
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


// ======================================================================
// conversion

template<typename F>
void cuda_mparticles_base::set_particles(uint n_prts, uint off, F getter, void *ctx)
{
  float4 *xi4  = new float4[n_prts];
  float4 *pxi4 = new float4[n_prts];
  
  for (int n = 0; n < n_prts; n++) {
    struct cuda_mparticles_prt prt = getter(n, ctx);

    for (int d = 0; d < 3; d++) {
      int bi = fint(prt.xi[d] * b_dxi[d]);
      if (bi < 0 || bi >= b_mx[d]) {
	printf("XXX xi %g %g %g\n", prt.xi[0], prt.xi[1], prt.xi[2]);
	printf("XXX n %d d %d xi4[n] %g biy %d // %d\n",
	       n, d, prt.xi[d], bi, b_mx[d]);
	if (bi < 0) {
	  prt.xi[d] = 0.f;
	} else {
	  prt.xi[d] *= (1. - 1e-6);
	}
      }
      bi = floorf(prt.xi[d] * b_dxi[d]);
      assert(bi >= 0 && bi < b_mx[d]);
    }

    xi4[n].x  = prt.xi[0];
    xi4[n].y  = prt.xi[1];
    xi4[n].z  = prt.xi[2];
    xi4[n].w  = cuda_int_as_float(prt.kind);
    pxi4[n].x = prt.pxi[0];
    pxi4[n].y = prt.pxi[1];
    pxi4[n].z = prt.pxi[2];
    pxi4[n].w = prt.qni_wni;
  }

  to_device(xi4, pxi4, n_prts, off);
  
  delete[] xi4;
  delete[] pxi4;
}

template<typename MP>
struct copy_ctx
{
  copy_ctx(MP& _mprts_other, int _p)
    : mprts_other(_mprts_other), p(_p)
  {
  }
  
  MP& mprts_other;
  int p;
};

template<typename MP>
struct ParticleGetter
{
  using particle_t = typename MP::particle_t;

  cuda_mparticles_prt operator()(int n, void *_ctx)
  {
    struct copy_ctx<MP> *ctx = (struct copy_ctx<MP> *) _ctx;
    const particle_t& prt_other = ctx->mprts_other[ctx->p][n];

    cuda_mparticles_prt prt;
    prt.xi[0]   = prt_other.xi;
    prt.xi[1]   = prt_other.yi;
    prt.xi[2]   = prt_other.zi;
    prt.pxi[0]  = prt_other.pxi;
    prt.pxi[1]  = prt_other.pyi;
    prt.pxi[2]  = prt_other.pzi;
    prt.kind    = prt_other.kind_;
    prt.qni_wni = prt_other.qni_wni;

    return prt;
  }
};

template<typename MP>
static void put_particle(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  using particle_t = typename MP::particle_t;
  struct copy_ctx<MP> *ctx = (struct copy_ctx<MP> *) _ctx;
  particle_t *part = &ctx->mprts_other[ctx->p][n];
  
  part->xi      = prt->xi[0];
  part->yi      = prt->xi[1];
  part->zi      = prt->xi[2];
  part->kind_   = prt->kind;
  part->pxi     = prt->pxi[0];
  part->pyi     = prt->pxi[1];
  part->pzi     = prt->pxi[2];
  part->qni_wni = prt->qni_wni;
}

template<typename MP>
static void copy_from(mparticles_cuda_t mprts, MP mprts_other)
{
  uint n_prts_by_patch[mprts->n_patches()];
  mprts->get_size_all(n_prts_by_patch);
  
  uint off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    copy_ctx<MP> ctx(mprts_other, p);
    ParticleGetter<MP> getter;
    mprts.sub_->cmprts()->set_particles(n_prts, off, getter, &ctx);

    off += n_prts;
  }
}

template<typename MP>
static void copy_to(mparticles_cuda_t mprts, MP mprts_other)
{
  uint n_prts_by_patch[mprts->n_patches()];
  mprts->get_size_all(n_prts_by_patch);
  
  uint off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    copy_ctx<MP> ctx(mprts_other, p);
    mprts->get_particles(n_prts, off, put_particle<MP>, &ctx);

    off += n_prts;
  }
}

// ======================================================================
// conversion to "single"

void psc_mparticles_cuda::copy_from_single(struct psc_mparticles *mprts_cuda,
					   struct psc_mparticles *mprts, uint flags)
{
  copy_from(mparticles_cuda_t(mprts_cuda), mparticles_single_t(mprts));
}

void psc_mparticles_cuda::copy_to_single(struct psc_mparticles *mprts_cuda,
					 struct psc_mparticles *mprts, uint flags)
{
  copy_to(mparticles_cuda_t(mprts_cuda), mparticles_single_t(mprts));
}

// ======================================================================
// conversion to "double"

void psc_mparticles_cuda::copy_from_double(struct psc_mparticles *mprts_cuda,
					   struct psc_mparticles *mprts, uint flags)
{
  copy_from(mparticles_cuda_t(mprts_cuda), mparticles_double_t(mprts));
}

void psc_mparticles_cuda::copy_to_double(struct psc_mparticles *mprts_cuda,
					 struct psc_mparticles *mprts, uint flags)
{
  copy_to(mparticles_cuda_t(mprts_cuda), mparticles_double_t(mprts));
}

