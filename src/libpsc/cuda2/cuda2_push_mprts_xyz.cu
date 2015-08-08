
#include "psc_cuda2.h"

#include "psc_particles_as_cuda2.h"
#include "psc_fields_cuda2.h"

#include "../psc_push_particles/inc_defs.h"

#define DIM DIM_XYZ
#define BLOCKSIZE_X 4
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4
#define EM_CACHE EM_CACHE_NONE
#define CALC_J CALC_J_1VB_VAR1

#define SFX(s) s ## _xyz

#include "cuda2_push_mprts.cu"

