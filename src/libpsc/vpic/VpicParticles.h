
#ifndef VPIC_PARTICLES_H
#define VPIC_PARTICLES_H

// ======================================================================
// VpicParticles

struct VpicSpecies : species_t
{
};

struct VpicParticles {
  species_t* append(species_t* s)
  {
    return ::append_species(s, reinterpret_cast<species_t **>(&sl_));
  }
  
  inline bool empty()
  {
    return !sl_;
  }

  
  //private:
  VpicSpecies *sl_;
};


#endif
