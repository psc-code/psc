
#ifndef VPIC_PARTICLES_BASE_H
#define VPIC_PARTICLES_BASE_H

#include "VpicListBase.h"

// ======================================================================
// VpicParticlesBase

struct VpicSpecies : species_t
{
};

struct VpicParticlesBase : public VpicListBase<VpicSpecies>
{
  typedef VpicSpecies Species;

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
  
  VpicSpecies* append(species_t* s)
  {
    species_t *sp = ::append_species(s, reinterpret_cast<species_t **>(&head_));
    return static_cast<VpicSpecies*>(sp);
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

  grid_t *getGrid_t()
  {
    assert(head_);
    return head_->g;
  }

  species_t* head()
  {
    return head_;
  }
};


#endif
