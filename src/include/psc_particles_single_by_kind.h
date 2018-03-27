
#ifndef PSC_PARTICLES_SINGLE_BY_KIND_H
#define PSC_PARTICLES_SINGLE_BY_KIND_H

#include "psc_particles_private.h"

#include "particles.hxx"

using particle_single_by_kind_real_t = float;

struct particle_single_by_kind_t
{
  using real_t = particle_single_by_kind_real_t;
};

struct MparticlesSingleByKind : MparticlesBase
{
  using Self = MparticlesSingleByKind;
  using particle_t = particle_single_by_kind_t;

  bk_mparticles *bkmprts;

  MparticlesSingleByKind(const Grid_t& grid)
    : MparticlesBase(grid)
  {
    bkmprts = bk_mparticles_new(grid.n_patches());
  }

  ~MparticlesSingleByKind()
  {
    bk_mparticles_delete(bkmprts);
  }

  static mrc_obj_method methods[];

  int get_n_prts() const override
  {
    return bk_mparticles_n_prts(bkmprts);
  }

  void get_size_all(uint *n_prts_by_patch) const override
  {
    bk_mparticles_size_all(bkmprts, n_prts_by_patch);
  }

  void reserve_all(const uint *n_prts_by_patch) override
  {
    bk_mparticles_reserve_all(bkmprts, n_prts_by_patch);
  }

  void resize_all(const uint *n_prts_by_patch) override
  {
    bk_mparticles_resize_all(bkmprts, n_prts_by_patch);
  }
};

// ======================================================================
// PscMparticlesSingleByKind

template<>
struct Mparticles_traits<MparticlesSingleByKind>
{
  static constexpr const char* name = "single_by_kind";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
