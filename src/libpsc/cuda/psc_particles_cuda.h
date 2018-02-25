
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "particles.hxx"
#include "particles_traits.hxx"
#include "psc_bits.h"

#include <vector>

#include "cuda_iface.h"

// ======================================================================
// PscMparticlesCuda

struct PscMparticlesCuda : PscMparticles<psc_mparticles_cuda>
{
  using Base = PscMparticles<psc_mparticles_cuda>;
  using particle_t = particle_cuda_t;
  using real_t = particle_cuda_real_t;
  using Real3 = Vec3<real_t>;
  using particle_buf_t = psc_particle_cuda_buf_t;

  using Base::Base;

  struct patch_t
  {
    patch_t(PscMparticlesCuda& mp, int p)
      : mp_(mp), p_(p), pi_(mp->grid())
    {
    }

    int blockPosition(real_t xi, int d) const { return pi_.blockPosition(xi, d); }
    Int3 blockPosition(const Real3& xi) const { return pi_.blockPosition(xi); }
    int validCellIndex(const particle_t& prt) const { return pi_.validCellIndex(&prt.xi); }
  
    const int* get_b_mx() const;

  private:
    PscMparticlesCuda& mp_;
    int p_;
    ParticleIndexer<real_t> pi_;
  };

  patch_t operator[](int p) {
    return patch_t(*this, p);
  }
};

template<>
struct mparticles_traits<PscMparticlesCuda>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};


#endif
