
#ifndef PSC_PARTICLES_BASE_H
#define PSC_PARTICLES_BASE_H

#include "VpicListBase.h"

inline void species_ctor(species_t *sp ,const char * name,
			 float q, float m, int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 grid_t *g)
{
}

// ======================================================================
// PscSpecies

template<class G>
struct PscSpecies : species_t
{
  typedef G Grid;
  
  // ----------------------------------------------------------------------
  // PscSpecies::ctor
  
  PscSpecies(const char * name, float q, float m,
	     int max_local_np, int max_local_nm,
	     int sort_interval, int sort_out_of_place,
	     Grid *grid)
  {
    CLEAR(this, 1);
  
    int len = name ? strlen(name) : 0;
    assert(len);
    assert(grid && grid->nv);
    assert(max_local_np > 0);
    assert(max_local_nm > 0);

    this->name = strdup(name);
    
    this->q = q;
    this->m = m;
    
    MALLOC_ALIGNED( this->p, max_local_np, 128 );
    this->max_np = max_local_np;
    
    MALLOC_ALIGNED( this->pm, max_local_nm, 128 );
    this->max_nm = max_local_nm;
    
    this->last_sorted       = INT64_MIN;
    this->sort_interval     = sort_interval;
    this->sort_out_of_place = sort_out_of_place;
    MALLOC_ALIGNED( this->partition, grid->nv+1, 128 );
    
    this->g = grid;
  }

  // ----------------------------------------------------------------------
  // PscSpecies::dtor
  
  ~PscSpecies()
  {
    FREE_ALIGNED(partition );
    FREE_ALIGNED(pm );
    FREE_ALIGNED(p );
    free(name);
  }
  
};

// ======================================================================
// PscParticlesBase

template<class G, class BCL>
struct PscParticlesBase : public VpicListBase<PscSpecies<G>>
{
  typedef G Grid;
  typedef BCL ParticleBcList;
  typedef PscSpecies<Grid> Species;
  typedef VpicListBase<Species> Base;

  using typename Base::iterator;
  using typename Base::const_iterator;
  using Base::size;
  using Base::begin;
  using Base::end;
  using Base::head_;
  using Base::push_front;

  static Species* create(const char * name, float q, float m,
			 int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 Grid *grid)
  {
    return new Species(name, q, m, max_local_np, max_local_nm,
		       sort_interval, sort_out_of_place, grid);
  }

  size_t getNumSpecies()
  {
    return size();
  }

  iterator find(int id)
  {
    return std::find_if(begin(), end(),
			[&id](const Species &sp) { return sp.id == id; });
  }
  
  iterator find(const char *name)
  {
    assert(name);
    return std::find_if(begin(), end(),
			[&name](const Species &sp) { return strcmp(sp.name, name) == 0; });
  }
  
  Species *append(Species *sp)
  {
    assert(!sp->next);
    if (find(sp->name) != end()) {
      ERROR(("There is already a material named \"%s\" in list", sp->name ));
    }
    int id = size();
    if (id >= ::max_material) {
      ERROR(("Too many materials in list to append material \"%s\"", sp->name));
    }
    sp->id   = (material_id)id;
    push_front(*sp);
    return sp;
  }
  
  grid_t *getGrid_t()
  {
    assert(head_);
    return head_->g;
  }

};
 
#endif
