
#ifndef PSC_SSE2_H
#define PSC_SSE2_H

////////
/// Toggle to switch precision 
///
/// 1 is double precision, 0 is single. 
#define SSE2_DOUBLE 1 

#include "psc.h"

// Kai's macros are much prettier (and a little faster) use them instead

#if 1

#define CF3(fldnr, jx,jy,jz)			\
  (sse2->flds.flds[(fldnr)*psc.fld_size + FF3_OFF(jx,jy,jz)])

#else
//out of range debugging
#define CF3(fldnr, jx,jy,jz)						\
  (*({int off = FF3_OFF(jx,jy,jz);					\
      assert(off >= 0);							\
      assert(off < psc.fld_size);					\
      &(sse2->flds.flds[(fldnr)*psc.fld_size + off]);			\
    }))

#endif


#include <assert.h>
#include <xmmintrin.h>
#include <emmintrin.h>
// Not including any SSE2 emulation at this time (finding an sse proc won't be hard, anything >= a P4 or AMD post 2005 will support these)

#include "simd_wrap.h"

/// Particle structure used by SSE2 implementation.
struct sse2_particle {
  sse2_real xi, yi, zi;
  sse2_real pxi, pyi, pzi;
  sse2_real qni;
  sse2_real mni;
  sse2_real wni;
};

/// Field/Particle data used by SSE2 Implementation
typedef struct psc_fields_sse2 {
  sse2_real *flds;
} psc_fields_sse2_t;

typedef struct psc_particles_sse2 {
  struct sse2_particle *particles; ///< Pointer to particle array
  int n_part;
} psc_particles_sse2_t;

struct psc_sse2 {
  psc_particles_sse2_t part;
  psc_fields_sse2_t flds;
  int part_allocated; ///< Number of particles currently allocated
};

// Packed vector datatypes, use typedefs to make things a bit prettier

/// Vector floating point type
typedef union packed_vector pvReal;

/// Vector integer type
typedef union packed_int pvInt;

/// VEC_SIZE particles packed into vectors
struct particle_vec{
  union packed_vector xi, yi, zi;
  union packed_vector pxi, pyi, pzi;
  union packed_vector qni, mni, wni;
};

void sse2_push_part_yz_a(void);
void sse2_push_part_yz_b(void);
void sse2_push_part_yz(void);
void init_vec_numbers(void);

// SSE2 needs to have these numbers packed into 
// vectors to utilize them effectively. 
pvReal ones, ///< Vector of "1.0"
  half, ///< Vector of "0.5"
  threefourths, ///< Vector of "0.75"
  onepfive, ///< Vector of "1.5"
  third; ///< Vector of "1./3."
pvInt ione; ///< Vector of "1"

void sse2_particles_from_fortran(struct psc_sse2 *sse2);
void sse2_particles_to_fortran(struct psc_sse2 *sse2);
void sse2_fields_from_fortran(struct psc_sse2 *sse2);
void sse2_fields_to_fortran(struct psc_sse2 *sse2);

#endif

//////////////////////////////////////////////////
/// \file psc_sse2.h SSE2 implementation datatypes and prototypes.
///
/// This file contains the macro SSE2_DOUBLE which sets the precision
/// (1 for double, 2 for single) of the sse2 code. The macros in 
/// simd_wrap.h choose the correct intrinsics based on this value. 
/// \FIXME It would be nice if precision in the C code could be set 
/// at compile time by some configure option. 
