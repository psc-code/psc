
#ifndef VPIC_PARTICLES_OPS
#define VPIC_PARTICLES_OPS

template<class P>
struct VpicParticlesOps
{
  typedef P Particles;
  typedef typename Particles::FieldArray FieldArray;
  typedef typename Particles::Interpolator Interpolator;
  typedef typename Particles::Accumulator Accumulator;
  
  VpicParticlesOps(vpic_simulation *simulation) : simulation_(simulation) { }

  void inject_particle(Particles& vmprts, Accumulator& accumulator, FieldArray& fa,
		       const struct psc_particle_inject *prt)
  {
    species_t *sp = &*vmprts.find(prt->kind);

    simulation_->inject_particle(sp, prt->x[0], prt->x[1], prt->x[2],
				 prt->u[0], prt->u[1], prt->u[2], prt->w, 0., 0);
  }

  void advance_p(Particles& vmprts, Accumulator& accumulator,
		 Interpolator& interpolator)
  {
    for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      TIC ::advance_p(&*sp, &accumulator, &interpolator); TOC(advance_p, 1);
    }
  }
  
  void boundary_p(const particle_bc_t& particle_bc_list, Particles& vmprts, FieldArray& fa,
		  Accumulator& accumulator)
  {
    ::boundary_p(const_cast<particle_bc_t*>(&particle_bc_list), vmprts.head(), &fa, &accumulator);
  }
  
  void accumulate_rho_p(Particles& vmprts, FieldArray &vmflds)
  {
    for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      TIC ::accumulate_rho_p(&vmflds, &*sp); TOC(accumulate_rho_p, 1);
    }
  }

  void accumulate_rhob(FieldArray& fa, const particle_t* p, float qsp)
  {
    ::accumulate_rhob(fa.f, p, fa.g, qsp);
  }

private:
  vpic_simulation *simulation_;
};


template<class ParticlesBase, class FA, class IA, class AA, class HA>
struct VpicParticles : ParticlesBase
{
  typedef ParticlesBase Base;
  typedef FA FieldArray;
  typedef IA Interpolator;
  typedef AA Accumulator;
  typedef HA HydroArray;

  using Base::Base;

  static void accumulate_hydro_p(HydroArray& ha, const species_t* sp,
				 const Interpolator& interpolator)
  {
    ::accumulate_hydro_p(&ha, sp, &interpolator);
  }

  static void uncenter_p(species_t *sp, const Interpolator& interpolator)
  {
    ::uncenter_p(sp, &interpolator);
  }
  
  static void sort_p(species_t *sp)
  {
    ::sort_p(sp);
  }

};

#endif

