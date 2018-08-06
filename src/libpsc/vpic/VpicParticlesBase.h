
#ifndef VPIC_PARTICLES_BASE_H
#define VPIC_PARTICLES_BASE_H

#include "VpicListBase.h"

// ======================================================================
// VpicParticlesBase

template<class G>
struct VpicSpecies : species_t
{
  typedef G Grid;
  
  Grid* grid() { return static_cast<Grid*>(g); }
  const Grid* grid() const { return static_cast<Grid*>(g); }
};

template<class G, class BCL>
struct VpicParticlesBase : public VpicListBase<VpicSpecies<G>>
{
  typedef G Grid;
  typedef BCL ParticleBcList;
  typedef VpicSpecies<Grid> Species;
  typedef VpicListBase<Species> Base;
  typedef particle_t Particle;
  typedef particle_mover_t ParticleMover;

  using typename Base::iterator;
  using typename Base::const_iterator;
  using Base::head_;

  static Species* create(const char * name, float q, float m,
			 int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 Grid *grid)
  {
    species_t* sp = species(name, q, m, max_local_np, max_local_nm,
			    sort_interval, sort_out_of_place, grid);
    return static_cast<Species*>(sp);
  }

  int getNumSpecies()
  {
    return ::num_species(head_);
  }
  
  Species* append(species_t* s)
  {
    species_t *sp = ::append_species(s, reinterpret_cast<species_t **>(&head_));
    return static_cast<Species*>(sp);
  }
  
  iterator find(int id)
  {
    species_t *sp = ::find_species_id(id, head_);
    return iterator(static_cast<Species*>(sp));
  }

  iterator find(const char *name)
  {
    species_t *sp = ::find_species_name(name, head_);
    return iterator(static_cast<Species*>(sp));
  }

  void inject_particle(VpicParticlesBase& vmprts, const particle_inject& prt)
  {
    species_t *sp = &*vmprts.find(prt.kind);

    extern vpic_simulation *simulation;
    assert(simulation);

    simulation->inject_particle(sp, prt.x[0], prt.x[1], prt.x[2],
				 prt.u[0], prt.u[1], prt.u[2], prt.w, 0., 0);
  }

  Grid *grid()
  {
    assert(head_);
    return head_->grid();
  }

  species_t* head()
  {
    return head_;
  }
};


#endif
