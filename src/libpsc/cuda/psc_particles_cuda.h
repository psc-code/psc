
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "particles.hxx"
#include "particles_traits.hxx"
#include "psc_bits.h"

#include <vector>

#include "cuda_iface.h"

// ======================================================================
// mparticles_base
// specialized for cuda to override sub_ handling (FIXME)

template<>
struct mparticles_base<psc_mparticles_cuda>
{
  using sub_t = psc_mparticles_cuda;
  using particle_t = typename sub_t::particle_t;
  using real_t = typename particle_t::real_t;

  explicit mparticles_base(psc_mparticles *mprts)
    : mprts_(mprts),
      sub_(mrc_to_subobj(mprts, sub_t))
  {
  }

  void put_as(psc_mparticles *mprts_base, unsigned int flags = 0)
  {
    psc_mparticles_put_as(mprts_, mprts_base, flags);
    mprts_ = nullptr;
  }

  unsigned int n_patches() { return mprts_->nr_patches; }
  
  psc_mparticles *mprts() { return mprts_; }
  
  sub_t* operator->() { return sub_; }

public:
  psc_mparticles *mprts_;
  sub_t *sub_;
};

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
};

template<>
struct mparticles_traits<mparticles_cuda_t>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};


#endif
