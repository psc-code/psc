
#ifndef SIMD_SSE2_H
#define SIMD_SSE2_H

#include <xmmintrin.h>

////////
/// Toggle to switch precision 
///
/// 1 is double precision, 0 is single. 

#define SSE2_DOUBLE 1

// ======================================================================

#if SSE2_DOUBLE

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 2

/// SSE2 floating point type
typedef double sse2_real;

/// SSE2 MPI type
#define MPI_SSE2_REAL MPI_DOUBLE

#else

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 4 

/// SSE2 floating point type
typedef float sse2_real;

/// SSE2 MPI type
#define MPI_SSE2_REAL MPI_FLOAT

#endif

#endif

