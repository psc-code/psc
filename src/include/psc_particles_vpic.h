
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  uint *n_prts_by_patch);

struct particle_vpic_t
{
  using real_t = float;
};

struct psc_mparticles_vpic : psc_mparticles_base
{
  using Base = psc_mparticles_base;
  using particle_t = particle_vpic_t; // FIXME, don't have it, but needed here...

  using Base::Base;
  
  Particles *vmprts;
  Simulation *sim;

  int get_n_prts() const override
  {
    return Simulation_mprts_get_nr_particles(sim, vmprts);
  }
  
  void get_size_all(uint *n_prts_by_patch) const override
  {
    vpic_mparticles_get_size_all(vmprts, 1, n_prts_by_patch);
  }

  void reserve_all(const uint *n_prts_by_patch) override
  {
    Simulation_mprts_reserve_all(sim, vmprts, 1, n_prts_by_patch);
  }

  void resize_all(const uint *n_prts_by_patch) override
  {
    Simulation_mprts_resize_all(sim, vmprts, 1, n_prts_by_patch);
  }

  void inject_reweight(int p, const psc_particle_inject *prt) override
  {
    Simulation_inject_particle(sim, vmprts, p, prt);
  }
};

using mparticles_vpic_t = mparticles_base<psc_mparticles_vpic>;

template<>
struct mparticles_traits<mparticles_vpic_t>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#define psc_mparticles_vpic(mprts) mrc_to_subobj(mprts, struct psc_mparticles_vpic)

#endif
