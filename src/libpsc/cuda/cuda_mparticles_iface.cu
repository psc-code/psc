
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

MparticlesCuda::MparticlesCuda(const Grid_t& grid)
  : MparticlesBase(grid)
{
  dprintf("CMPRTS: ctor\n");
  cmprts_ = new cuda_mparticles(grid);

  patches_.reserve(grid.n_patches());
  for (int p = 0; p < grid.n_patches(); p++) {
      patches_.emplace_back(*this, p);
  }
}

MparticlesCuda::~MparticlesCuda()
{
  dprintf("CMPRTS: dtor\n");
  delete cmprts_;
}

void MparticlesCuda::reserve_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: reserve_all\n");
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  cmprts_->reserve_all(n_prts_by_patch);
}

void MparticlesCuda::get_size_all(uint *n_prts_by_patch) const
{
  dprintf("CMPRTS: get_size_all\n");
  cmprts_->get_size_all(n_prts_by_patch);
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
}

void MparticlesCuda::resize_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: resize_all\n");
  cmprts_->resize_all(n_prts_by_patch);
}

int MparticlesCuda::get_n_prts() const
{
  dprintf("CMPRTS: get_n_prts\n");
  return cmprts_->get_n_prts();
}

void MparticlesCuda::to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off)
{
  dprintf("CMPRTS: to_device\n");
  cmprts_->to_device(xi4, pxi4, n_prts, off);
}

void MparticlesCuda::from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off)
{
  dprintf("CMPRTS: from_device\n");
  cmprts_->from_device(xi4, pxi4, n_prts, off);
}

void MparticlesCuda::setup_internals()
{
  dprintf("CMPRTS: setup_internals\n");
  cmprts_->setup_internals();
}

void MparticlesCuda::inject_buf(cuda_mparticles_prt *buf, uint *buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
  cmprts_->inject(buf, buf_n_by_patch);
}

const int* MparticlesCuda::patch_get_b_mx(int p)
{
  return cmprts_->patch_get_b_mx(p);
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
    prt.qni_wni = prt_other.qni_wni_;

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
    prt_other.qni_wni_ = prt.qni_wni;
  }

private:
  MP& mprts_other_;
  int p_;
};

template<typename MP>
void MparticlesCuda::copy_from(MparticlesCuda& mp, MP& mp_other)
{
  int n_patches = mp.n_patches();
  uint n_prts_by_patch[n_patches];
  mp_other.get_size_all(n_prts_by_patch);
  mp.reserve_all(n_prts_by_patch);
  mp.resize_all(n_prts_by_patch);

  for (int p = 0; p < n_patches; p++) {
    ConvertToCuda<MP> convert_to_cuda(mp_other, p);
    mp.cmprts()->set_particles(p, convert_to_cuda);
  }

  mp.setup_internals();
}

template<typename MP>
void MparticlesCuda::copy_to(MparticlesCuda& mp, MP& mp_other)
{
  int n_patches = mp_other.n_patches();
  uint n_prts_by_patch[n_patches];
  mp.get_size_all(n_prts_by_patch);
  mp_other.reserve_all(n_prts_by_patch);
  mp_other.resize_all(n_prts_by_patch);

  for (int p = 0; p < n_patches; p++) {
    ConvertFromCuda<MP> convert_from_cuda(mp_other, p);
    mp.cmprts()->get_particles(p, convert_from_cuda);
  }
}

// ======================================================================
// conversion to "single"/"double"

template<typename MP>
void MparticlesCuda::copy_from(MparticlesBase& mp, MparticlesBase& mp_other)
{
  copy_from(dynamic_cast<MparticlesCuda&>(mp), dynamic_cast<MP&>(mp_other));
}

template<typename MP>
void MparticlesCuda::copy_to(MparticlesBase& mp, MparticlesBase& mp_other)
{
  copy_to(dynamic_cast<MparticlesCuda&>(mp), dynamic_cast<MP&>(mp_other));
}

const MparticlesCuda::Convert MparticlesCuda::convert_to_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_to<MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_to<MparticlesDouble>   },
};

const MparticlesCuda::Convert MparticlesCuda::convert_from_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_from<MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_from<MparticlesDouble>   },
};

