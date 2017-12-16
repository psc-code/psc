
#ifndef PSC_PARTICLES_BASE_H
#define PSC_PARTICLES_BASE_H

#include "VpicListBase.h"

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
    
    this->p = new particle_t[max_local_np];
    this->max_np = max_local_np;
    
    this->pm = new particle_mover_t[max_local_nm];
    this->max_nm = max_local_nm;
    
    this->last_sorted       = INT64_MIN;
    this->sort_interval     = sort_interval;
    this->sort_out_of_place = sort_out_of_place;
    this->partition = new int[grid->nv + 1];
    
    this->g = grid;
  }

  // ----------------------------------------------------------------------
  // PscSpecies::dtor
  
  ~PscSpecies()
  {
    delete[] partition;
    delete[] pm;
    delete[] p;
    free(name);
  }

  Grid* getGrid()
  {
    return static_cast<Grid*>(g);
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
  
  Grid *getGrid()
  {
    assert(head_);
    return reinterpret_cast<Grid*>(head_->g);
  }
  
  Species* head()
  {
    return head_;
  }
};
 
#endif
