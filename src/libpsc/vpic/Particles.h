
#ifndef PARTICLES_H
#define PARTICLES_H

// ======================================================================
// Particles

struct Particles {
  Particles(species_t*& sl);

  species_t* append(species_t* sl);
  bool empty();
  
  //private:
  species_t *&sl_;
};

inline Particles::Particles(species_t*& sl)
  : sl_(sl)
{
}

inline species_t* Particles::append(species_t* s)
{
  return ::append_species(s, &sl_);
}

inline bool Particles::empty()
{
  return !sl_;
}

#endif
