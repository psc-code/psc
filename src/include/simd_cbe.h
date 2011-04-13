
#ifndef SIMD_CBE_H
#define SIMD_CBE_H

#define CBE_DOUBLE 1

// ======================================================================

#if CBE_DOUBLE

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 2

/// CBE floating point type
typedef double cbe_real;

/// CBE MPI type
#define MPI_CBE_REAL MPI_DOUBLE

#else

/// Number of elements in a floating point vector (changes with precision)
#define VEC_SIZE 4 

/// CBE floating point type
typedef float cbe_real;

/// CBE MPI type
#define MPI_CBE_REAL MPI_FLOAT

#endif

#endif

