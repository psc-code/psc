
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "particles.hxx"
#include "particles_traits.hxx"

#include "psc_particle_buf_cuda.h"

// ======================================================================
// psc_mparticles_cuda

struct psc_mparticles_cuda {
  struct cuda_mparticles *cmprts;
};

#define psc_mparticles_cuda(mprts) mrc_to_subobj(mprts, struct psc_mparticles_cuda)

// ======================================================================
// mparticles_cuda_t

struct mparticles_cuda_t : mparticles_base<psc_mparticles_cuda>
{
  using Base = mparticles_base<psc_mparticles_cuda>;
  using particle_t = particle_cuda_t;
  using real_t = particle_cuda_real_t;
  using particle_buf_t = psc_particle_cuda_buf_t;

  using Base::Base;
  
  struct patch_t
  {
    patch_t(mparticles_cuda_t& mp, int p)
      : mp_(mp), p_(p)
    {
    }

    const int* get_b_mx() const;
    const real_t* get_b_dxi() const;
    
  private:
    mparticles_cuda_t& mp_;
    int p_;
  };

  patch_t operator[](int p) {
    return patch_t(*this, p);
  }
};

template<>
struct mparticles_traits<mparticles_cuda_t>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

template<>
struct mparticles_traits<particle_cuda_t>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};



#endif
