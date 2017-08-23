
#ifndef PSC_SSE2_H
#define PSC_SSE2_H

#include "psc.h"
#include "psc_particles_sse2.h"
#include "psc_fields_sse2.h"

#include <assert.h>
#include <emmintrin.h>
// Not including any SSE2 emulation at this time (finding an sse proc won't be hard, anything >= a P4 or AMD post 2005 will support these)

#include "simd_wrap.h"

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

void psc_push_particles_sse2_push_yz(struct psc_push_particles *push,
				     struct psc_mparticles *particles_base,
				     struct psc_mfields *flds_base);
void psc_push_particles_sse2_push_yz_a(struct psc_push_particles *push,
				       struct psc_mparticles *particles_base,
				       struct psc_mfields *flds_base);
void psc_push_particles_sse2_push_yz_b(struct psc_push_particles *push,
				       struct psc_mparticles *particles_base,
				       struct psc_mfields *flds_base);
void init_vec_numbers(void);
__m128i func_mul_epu32(__m128i a, __m128i b);

// SSE2 needs to have these numbers packed into 
// vectors to utilize them effectively. 
pvReal ones, ///< Vector of "1.0"
  half, ///< Vector of "0.5"
  threefourths, ///< Vector of "0.75"
  onepfive, ///< Vector of "1.5"
  third; ///< Vector of "1./3."
pvInt ione; ///< Vector of "1"

#endif

//////////////////////////////////////////////////
/// \file psc_sse2.h SSE2 implementation datatypes and prototypes.
///
/// This file contains the macro SSE2_DOUBLE which sets the precision
/// (1 for double, 2 for single) of the sse2 code. The macros in 
/// simd_wrap.h choose the correct intrinsics based on this value. 
/// \FIXME It would be nice if precision in the C code could be set 
/// at compile time by some configure option. 
