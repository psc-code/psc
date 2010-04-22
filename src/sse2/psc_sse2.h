
#ifndef PSC_SSE2_H
#define PSC_SSE2_H

#include <assert.h>
#include <xmmintrin.h>
// Not including any SSE2 emulation at this time (finding an sse proc won't be hard, anything >= a P4 or AMD post 2005 will support these)

#include "psc.h"
// For now, we'll stick with the Array of Structs memory layout, and do our packing at the loop level.
// There is some question as to whether packing when the particles are brought in will have any impact on performance...

struct sse2_particle {
  real xi, yi, zi;
  real pxi, pyi, pzi;
  real qni;
  real mni;
  real wni;
};

struct psc_sse2 {
  struct sse2_particle *part;
};

void sse2_push_part_yz_a();
void sse2_particles_from_fortran(struct psc_sse2 *sse2);
void sse2_particles_to_fortran(struct psc_sse2 *sse2);

#endif
