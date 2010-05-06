
#ifndef PSC_SSE2_H
#define PSC_SSE2_H

#define SSE2_FLOAT  // (or SSE2_DOUBLE) change this to switch between float and double sse2 

#define FORT_FIELD(T,l,j,k) (float)psc.f_fields[(T)][((l)-psc.ilo[0]+psc.ibn[0]) + ((j)-psc.ilo[1]+psc.ibn[1])*(psc.img[0]) + ((k) - psc.ilo[2]+psc.ibn[2])*(psc.img[0]*psc.img[1])]

#define C_FIELD(T,l,j,k) sse2->fields[(T)*psc.fld_size + ((l)-psc.ilo[0]+psc.ibn[0]) + ((j)-psc.ilo[1]+psc.ibn[1])*(psc.img[0]) + ((k) - psc.ilo[2]+psc.ibn[2])*(psc.img[0]*psc.img[1])]

#include <assert.h>
#include <xmmintrin.h>
#include <emmintrin.h>
// Not including any SSE2 emulation at this time (finding an sse proc won't be hard, anything >= a P4 or AMD post 2005 will support these)

#include "psc.h"

#include "sse2_cgen.h"

struct sse2_particle {
  sse2_real xi, yi, zi;
  sse2_real pxi, pyi, pzi;
  sse2_real qni;
  sse2_real mni;
  sse2_real wni;
};

struct psc_sse2 {
  struct sse2_particle *part;
  sse2_real *fields;
};

// Packed vector datatypes, use typedefs to make things a bit prettier

typedef union packed_vector pvReal;

typedef union packed_int pvInt;


void sse2_push_part_yz_a();
void sse2_push_part_yz_b();
void sse2_push_part_yz();

#endif
