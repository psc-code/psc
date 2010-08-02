
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc.h"

// this matches the Fortran particle data structure

struct f_particle {
  f_real xi, yi, zi;
  f_real pxi, pyi, pzi;
  f_real qni;
  f_real mni;
  f_real cni;
  f_real lni;
  f_real wni;
};



#endif
