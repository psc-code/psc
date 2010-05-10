
#ifndef PSC_SSE2_H
#define PSC_SSE2_H

#define SSE2_DOUBLE 0 // toggle to switch precision 
#include "psc.h"

// Kai's macros are much prettier (and a little faster) use them instead

#if 1

#define CF3(fldnr, jx,jy,jz)			\
  (sse2->fields[(fldnr)*psc.fld_size + FF3_OFF(jx,jy,jz)])

#else
//out of range debugging
#define CF3(fldnr, jx,jy,jz)						\
  (*({int off = FF3_OFF(jx,jy,jz);					\
      assert(off >= 0);							\
      assert(off < psc.fld_size);					\
      &(sse2->fields[(fldnr)*psc.fld_size + off]);					\
    }))

#endif


#include <assert.h>
#include <xmmintrin.h>
#include <emmintrin.h>
// Not including any SSE2 emulation at this time (finding an sse proc won't be hard, anything >= a P4 or AMD post 2005 will support these)

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

struct particle_vec{
  union packed_vector xi, yi, zi;
  union packed_vector pxi, pyi, pzi;
  union packed_vector qni, mni, wni;
};

void sse2_push_part_yz_a();
void sse2_push_part_yz_b();
void sse2_push_part_yz();

#endif
