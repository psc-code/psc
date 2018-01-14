
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "particles.hxx"
#include "particles_traits.hxx"

#include "psc_particle_buf_cuda.h"

// ======================================================================
// psc_mparticles_cuda

// struct psc_mparticles_cuda {
//   using particle_t = particle_cuda_t;
//   struct cuda_mparticles *cmprts;
// };

// #define psc_mparticles_cuda(mprts) mrc_to_subobj(mprts, struct psc_mparticles_cuda)

// ======================================================================
// mparticles_base

// template<>
// struct mparticles_base<psc_mparticles_cuda>
// {
//   using sub_t = psc_mparticles_cuda;
//   using particle_t = particle_cuda_t;
//   using real_t = typename particle_t::real_t;

//   explicit mparticles_base(psc_mparticles *mprts) : mprts_(mprts) { }

//   void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
//   {
//     psc_mparticles_put_as(mprts_, mprts_base, flags);
//   }

//   unsigned int n_patches() { return mprts_->nr_patches; }
  
//   psc_mparticles *mprts() { return mprts_; }
  
// private:
//   psc_mparticles *mprts_;
// };

// ======================================================================
// mparticles_cuda_t

struct mparticles_cuda_t
{
  using particle_t = particle_cuda_t;
  using real_t = particle_cuda_real_t;
  using particle_buf_t = psc_particle_cuda_buf_t;

  struct cuda_mparticles* cmprts()
  {
    return cmprts_;
  }
  
  struct patch_t
  {
    patch_t(mparticles_cuda_t& mp, int p)
      : mp_(mp), p_(p)
    {
    }

    void get_block_pos(const real_t xi[3], int b_pos[3])
    {
      const real_t* b_dxi = get_b_dxi();
      for (int d = 0; d < 3; d++) {
	b_pos[d] = fint(xi[d] * b_dxi[d]);
      }
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

  mparticles_cuda_t* sub() { return this; }

  void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
  {
    psc_mparticles_put_as(mprts_, mprts_base, flags);
  }

  unsigned int n_patches() { return mprts_->nr_patches; }
  
  psc_mparticles *mprts() { return mprts_; }
  
  // private:
  psc_mparticles *mprts_;
  cuda_mparticles* cmprts_;
};

using psc_mparticles_cuda = mparticles_cuda_t;

inline mparticles_cuda_t& _mparticles_cuda_t(psc_mparticles* mprts)
{
  auto rv = mrc_to_subobj(mprts, struct mparticles_cuda_t);
  rv->mprts_ = mprts;
  return *rv;
}

template<>
struct mparticles_traits<mparticles_cuda_t>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

template<>
inline mparticles_cuda_t& psc_mparticles::get_as(unsigned int flags)
{
  const char *type = mparticles_traits<mparticles_cuda_t>::name;
  struct psc_mparticles *mprts = psc_mparticles_get_as(this, type, flags);
  return _mparticles_cuda_t(mprts);
}


#endif
