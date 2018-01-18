
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

template<typename MP>
struct ConvertToCuda
{
  using particle_t = typename MP::particle_t;

  ConvertToCuda(MP& mprts_other, int p)
    : mprts_other_(mprts_other), p_(p)
  {
  }

  cuda_mparticles_prt operator()(int n)
  {
    const particle_t& prt_other = mprts_other_[p_][n];

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

private:
  MP& mprts_other_;
  int p_;
};

template<typename MP>
struct ConvertFromCuda
{
  using particle_t = typename MP::particle_t;

  ConvertFromCuda(MP& mprts_other, int p)
    : mprts_other_(mprts_other), p_(p)
  {
  }

  void operator()(int n, const cuda_mparticles_prt &prt)
  {
    particle_t& prt_other = mprts_other_[p_][n];

    prt_other.xi      = prt.xi[0];
    prt_other.yi      = prt.xi[1];
    prt_other.zi      = prt.xi[2];
    prt_other.kind_   = prt.kind;
    prt_other.pxi     = prt.pxi[0];
    prt_other.pyi     = prt.pxi[1];
    prt_other.pzi     = prt.pxi[2];
    prt_other.qni_wni = prt.qni_wni;
  }

private:
  MP& mprts_other_;
  int p_;
};

template<typename MP>
static void copy_from(mparticles_cuda_t mprts_to, MP mprts_from)
{
  int n_patches = mprts_to.n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from.get_size_all((int*) n_prts_by_patch); // FIXME cast
  mprts_to.reserve_all((int*) n_prts_by_patch); // FIXME cast
  mprts_to.resize_all((int*) n_prts_by_patch); // FIXME cast

  mprts_to->get_size_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertToCuda<MP> convert_to_cuda(mprts_from, p);
    mprts_to->cmprts()->set_particles(n_prts, off, convert_to_cuda);

    off += n_prts;
  }
}

template<typename MP>
static void copy_to(mparticles_cuda_t mprts_from, MP mprts_to)
{
  int n_patches = mprts_to.n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from.get_size_all((int*) n_prts_by_patch); // FIXME cast
  mprts_to.reserve_all((int*) n_prts_by_patch); // FIXME cast
  mprts_to.resize_all((int*) n_prts_by_patch); // FIXME cast

  mprts_from->get_size_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    ConvertFromCuda<MP> convert_from_cuda(mprts_to, p);
    mprts_from->cmprts()->get_particles(n_prts, off, convert_from_cuda);

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

