
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include "psc_particle_buf_cuda.h"

#define PTYPE PTYPE_CUDA
#include "psc_particles_common.h"
#undef PTYPE

// ======================================================================
// psc_mparticles_cuda

struct psc_mparticles_cuda {
  struct cuda_mparticles *cmprts;
};

const particle_cuda_real_t *psc_mparticles_cuda_patch_get_b_dxi(struct psc_mparticles *mprts, int p);
const int *psc_mparticles_cuda_patch_get_b_mx(struct psc_mparticles *mprts, int p);

// ======================================================================
// mparticles_cuda_t

struct mparticles_cuda_t : mparticles_base<psc_mparticles_cuda>
{
  using Base = mparticles_base<psc_mparticles_cuda>;
  using real_t = particle_cuda_real_t;

  using Base::Base;
  
  struct patch_t
  {
    patch_t(mparticles_cuda_t& mp, int p)
      : mp_(mp), p_(p)
    {
    }

    const int* get_b_mx() const
    {
      return psc_mparticles_cuda_patch_get_b_mx(mp_.mprts(), p_);
    }
    
    const real_t* get_b_dxi() const
    {
      return psc_mparticles_cuda_patch_get_b_dxi(mp_.mprts(), p_);
    }
    
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
