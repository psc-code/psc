
#ifndef VPIC_PARTICLES_BASE_H
#define VPIC_PARTICLES_BASE_H

#include "VpicListBase.h"

// ======================================================================
// VpicParticlesBase

inline void species_ctor(species_t *sp ,const char * name,
			 float q, float m, int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 grid_t *grid);
inline void species_dtor(species_t *sp);

struct VpicSpecies : species_t
{
  VpicSpecies(const char * name, float q, float m,
	      int max_local_np, int max_local_nm,
	      int sort_interval, int sort_out_of_place,
	      Grid *grid)
  {
    species_ctor(this, name, q, m, max_local_np, max_local_nm,
		 sort_interval, sort_out_of_place, grid);
  }

  ~VpicSpecies()
  {
    species_dtor(this);
  }
  
};

struct VpicParticlesBase : public VpicListBase<VpicSpecies>
{
  typedef VpicSpecies Species;

  int getNumSpecies()
  {
    return ::num_species(head_);
  }
  
  VpicSpecies* append(species_t* s)
  {
    return static_cast<VpicSpecies*>(::append_species(s, reinterpret_cast<species_t **>(&head_)));
  }
  
  bool empty()
  {
    return !head_;
  }

  int size()
  {
    int sz = 0;
    for (const_iterator sp = cbegin(); sp != cend(); ++sp) {
      sz++;
    }
    return sz;
  }

  iterator find_id(int id)
  {
    for (auto sp = begin(); sp != end(); ++sp) {
      if (sp->id == id) {
	return sp;
      }
    }
    return end();
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

// ======================================================================
// c&p

inline void species_ctor(species_t *sp ,const char * name,
			 float q, float m, int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 grid_t *g)
{
  CLEAR(sp, 1);
  
  int len = name ? strlen(name) : 0;

  if( !len ) ERROR(( "Cannot create a nameless species" ));
  if( !g ) ERROR(( "NULL grid" ));
  if( g->nv == 0) ERROR(( "Allocate grid before defining species." ));
  if( max_local_np<1 ) max_local_np = 1;
  if( max_local_nm<1 ) max_local_nm = 1;

  MALLOC( sp->name, len+1 );
  strcpy( sp->name, name );

  sp->q = q;
  sp->m = m;

  MALLOC_ALIGNED( sp->p, max_local_np, 128 );
  sp->max_np = max_local_np;

  MALLOC_ALIGNED( sp->pm, max_local_nm, 128 );
  sp->max_nm = max_local_nm;

  sp->last_sorted       = INT64_MIN;
  sp->sort_interval     = sort_interval;
  sp->sort_out_of_place = sort_out_of_place;
  MALLOC_ALIGNED( sp->partition, g->nv+1, 128 );

  sp->g = g;
}

inline void species_dtor(species_t *sp)
{
  FREE_ALIGNED( sp->partition );
  FREE_ALIGNED( sp->pm );
  FREE_ALIGNED( sp->p );
  FREE( sp->name );
}

#endif
