
#include "psc.h"

// ----------------------------------------------------------------------
// generic C data structures

struct c_particle {
  real xi, yi, zi;
  real pxi, pyi, pzi;
  real qni;
  real mni;
  real wni;
};

struct psc_genc {
  struct c_particle *part;
};

void genc_push_part_yz_a();


