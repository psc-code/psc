
#define PTYPE_SINGLE          1
#define PTYPE_DOUBLE          2
#define PTYPE_SINGLE_BY_BLOCK 3
#define PTYPE_C               4
#define PTYPE_FORTRAN         5

#if PTYPE == PTYPE_SINGLE

#define particle_real_t particle_single_real_t

#elif PTYPE == PTYPE_DOUBLE

#define particle_real_t particle_double_real_t

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define particle_real_t particle_single_by_block_real_t

#elif PTYPE == PTYPE_C

#define particle_real_t particle_c_real_t

#elif PTYPE == PTYPE_FORTRAN

#define particle_real_t particle_fortran_real_t

#endif


#if PTYPE == PTYPE_SINGLE

typedef float particle_real_t;
#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

#elif PTYPE == PTYPE_DOUBLE

typedef double particle_real_t;
#define MPI_PARTICLES_DOUBLE_REAL MPI_DOUBLE

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

typedef float particle_real_t;
#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

#elif PTYPE == PTYPE_C

typedef double particle_real_t;
#define MPI_PARTICLES_C_REAL MPI_DOUBLE

#elif PTYPE == PTYPE_FORTRAN

typedef double particle_real_t;
#define MPI_PARTICLES_FORTRAN_REAL MPI_DOUBLE

#endif

#include <math.h>

#undef particle_real_t

