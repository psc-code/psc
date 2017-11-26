
#ifndef VPIC_PARTICLES_OPS
#define VPIC_PARTICLES_OPS

template<class P, class IA, class AA, class FA>
struct VpicParticlesOps {
  typedef P Particles;
  typedef FA FieldArray;
  typedef IA Interpolator;
  typedef AA Accumulator;
  
  VpicParticlesOps(vpic_simulation *simulation) : simulation_(simulation) { }

  void inject_particle(Particles& vmprts, Accumulator& accumulator, FieldArray& fa,
		       const struct psc_particle_inject *prt)
  {
    species_t *sp = &*vmprts.find_id(prt->kind);

    simulation_->inject_particle(sp, prt->x[0], prt->x[1], prt->x[2],
				 prt->u[0], prt->u[1], prt->u[2], prt->w, 0., 0);
  }

  void advance_p(Particles& vmprts, Accumulator& accumulator,
		 Interpolator& interpolator)
  {
    for (typename Particles::Iter sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      TIC ::advance_p(&*sp, &accumulator, &interpolator); TOC(advance_p, 1);
    }
  }
  
  void boundary_p(particle_bc_t *particle_bc_list, Particles& vmprts, FieldArray& fa,
		  Accumulator& accumulator)
  {
    ::boundary_p(particle_bc_list, vmprts.head(), &fa, &accumulator);
  }
  
  void accumulate_rho_p(Particles& vmprts, FieldArray &vmflds)
  {
    for (typename Particles::Iter sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
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

#endif

