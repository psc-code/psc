
#ifndef VPIC_PARTICLES_BASE_H
#define VPIC_PARTICLES_BASE_H

// ======================================================================
// VpicParticlesBase

struct VpicSpecies : species_t
{
};

struct VpicSpeciesIter {
  typedef VpicSpeciesIter Iter;
  
  VpicSpeciesIter(VpicSpecies *node=0) : node_(node)
  {
  }

  bool operator!=(const Iter& x) const
  {
    return node_ != x.node_;
  }

  Iter& operator++()
  {
    node_ = static_cast<VpicSpecies *>(node_->next);
    return *this;
  }

  VpicSpecies& operator*() const
  {
    return *node_;
  }
  
  VpicSpecies *operator->() const
  {
    return node_;
  }
  
private:
  VpicSpecies *node_;
};

struct VpicParticlesBase {
  typedef VpicSpecies Species;
  typedef VpicSpeciesIter Iter;

  int getNumSpecies()
  {
    return ::num_species(sl_);
  }
  
  species_t* append(species_t* s)
  {
    return ::append_species(s, reinterpret_cast<species_t **>(&sl_));
  }
  
  bool empty()
  {
    return !sl_;
  }

  int size()
  {
    int sz = 0;
    for (Iter sp = begin(); sp != end(); ++sp) {
      sz++;
    }
    return sz;
  }

  VpicSpeciesIter begin()
  {
    return VpicSpeciesIter(sl_);
  }
  
  VpicSpeciesIter end()
  {
    return VpicSpeciesIter();
  }

  VpicSpeciesIter find_id(int id)
  {
    Iter sp;
    for (sp = begin(); sp != end(); ++sp) {
      if (sp->id == id) {
	break;
      }
    }
    return sp;
  }

  grid_t *getGrid_t()
  {
    return sl_->g;
  }

  species_t* head()
  {
    return sl_;
  }
  
private:
  VpicSpecies *sl_;
};


#endif
