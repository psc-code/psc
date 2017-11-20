
#ifndef PARTICLES_H
#define PARTICLES_H

// ======================================================================
// VpicParticles

struct VpicParticles {
  VpicParticles(species_t*& sl);

  species_t* append(species_t* sl);
  bool empty();
  
  //private:
  species_t *&sl_;
};

inline VpicParticles::VpicParticles(species_t*& sl)
  : sl_(sl)
{
}

inline species_t* VpicParticles::append(species_t* s)
{
  return ::append_species(s, &sl_);
}

inline bool VpicParticles::empty()
{
  return !sl_;
}

#endif
