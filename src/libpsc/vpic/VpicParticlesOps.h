
#ifndef VPIC_PARTICLES_OPS
#define VPIC_PARTICLES_OPS

template<class P>
struct VpicParticlesOps {
  typedef P Particles;
  
  VpicParticlesOps(vpic_simulation *simulation) : simulation_(simulation) { }

  void inject_particle(Particles *vmprts, int p, const struct psc_particle_inject *prt)
  {
    assert(p == 0);
    species_t *sp = &*vmprts->find_id(prt->kind);

    simulation_->inject_particle(sp, prt->x[0], prt->x[1], prt->x[2],
				 prt->u[0], prt->u[1], prt->u[2], prt->w, 0., 0);
  }

private:
  vpic_simulation *simulation_;
};

#endif

