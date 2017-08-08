
#define PTYPE_SINGLE          1
#define PTYPE_DOUBLE          2
#define PTYPE_SINGLE_BY_BLOCK 3
#define PTYPE_C               4
#define PTYPE_FORTRAN         5

#if PTYPE == PTYPE_SINGLE

#define particle_real_t particle_single_real_t
#define particle_t particle_single_t
#define psc_particle psc_particle_single

#elif PTYPE == PTYPE_DOUBLE

#define particle_real_t particle_double_real_t
#define particle_t particle_double_t
#define psc_particle psc_particle_double

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define particle_real_t particle_single_by_block_real_t
#define particle_t particle_single_by_block_t
#define psc_particle psc_particle_single_by_block

#elif PTYPE == PTYPE_C

#define particle_real_t particle_c_real_t
#define particle_t particle_c_t
#define psc_particle psc_particle_c

#elif PTYPE == PTYPE_FORTRAN

#define particle_real_t particle_fortran_real_t
#define particle_t particle_fortran_t
#define psc_particle psc_particle_fortran

#endif

// ----------------------------------------------------------------------
// particle_real_t

#if PTYPE == PTYPE_SINGLE || PTYPE == PTYPE_SINGLE_BY_BLOCK

typedef float particle_real_t;

#elif PTYPE == PTYPE_DOUBLE || PTYPE == PTYPE_C || PTYPE == PTYPE_FORTRAN

typedef double particle_real_t;

#endif

// ----------------------------------------------------------------------
// MPI_PARTICLES_REAL
// annoying, but need to use a macro, which means we can't consolidate float/double

#if PTYPE == PTYPE_SINGLE

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

#elif PTYPE == PTYPE_DOUBLE

#define MPI_PARTICLES_DOUBLE_REAL MPI_DOUBLE

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

#elif PTYPE == PTYPE_C

#define MPI_PARTICLES_C_REAL MPI_DOUBLE

#elif PTYPE == PTYPE_FORTRAN

#define MPI_PARTICLES_FORTRAN_REAL MPI_DOUBLE

#endif

// ----------------------------------------------------------------------
// particle_t

#if PTYPE == PTYPE_SINGLE || PTYPE == PTYPE_SINGLE_BY_BLOCK || PTYPE == PTYPE_DOUBLE

typedef struct psc_particle {
  particle_real_t xi, yi, zi;
  particle_real_t qni_wni;
  particle_real_t pxi, pyi, pzi;
  int kind;
} particle_t;

#elif PTYPE == PTYPE_C

typedef struct psc_particle {
  particle_real_t xi, yi, zi;
  particle_real_t pxi, pyi, pzi;
  particle_real_t qni;
  particle_real_t mni;
  particle_real_t wni;
  long long kind; // 64 bits to match the other members, for bnd exchange
} particle_t;

#elif PTYPE == PTYPE_FORTRAN

typedef struct psc_particle {
  particle_real_t xi, yi, zi;
  particle_real_t pxi, pyi, pzi;
  particle_real_t qni;
  particle_real_t mni;
  particle_real_t cni;
  particle_real_t lni;
  particle_real_t wni;
} particle_t;

#endif

// ----------------------------------------------------------------------
// psc_mparticles_patch

#if PTYPE == PTYPE_SINGLE

struct psc_mparticles_single_patch {
  particle_single_t *prt_array;
  int n_prts;
  int n_alloced;

  particle_single_t *prt_array_alt;
  int b_mx[3];
  int nr_blocks;
  particle_single_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int n_send;
  unsigned int n_part_save;
  bool need_reorder;
};

#elif PTYPE == PTYPE_DOUBLE

struct psc_mparticles_double_patch {
  particle_double_t *prt_array;
  int n_prts;
  int n_alloced;
};

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

struct psc_mparticles_single_by_block_patch {
  particle_single_by_block_t *prt_array;
  int n_prts;
  int n_alloced;

  particle_single_by_block_t *prt_array_alt;
  int b_mx[3];
  int nr_blocks;
  particle_single_by_block_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;
  bool need_reorder;
};

#elif PTYPE == PTYPE_C

struct psc_mparticles_c_patch {
  particle_c_t *prt_array;
  int n_prts;
  int n_alloced;
};

#elif PTYPE == PTYPE_FORTRAN

struct psc_mparticles_fortran_patch {
  particle_fortran_t *prt_array;
  int n_prts;
  int n_alloced;
};

#endif

// ----------------------------------------------------------------------
//

#if PTYPE == PTYPE_SINGLE

#elif PTYPE == PTYPE_DOUBLE

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#elif PTYPE == PTYPE_C

#elif PTYPE == PTYPE_FORTRAN

#endif




#include <math.h>

#undef particle_real_t
#undef particle_t
#undef psc_particle
