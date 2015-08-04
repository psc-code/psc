
#include "psc_debug.h"

#include "psc_particles_as_cuda2.h"

#include "psc_fields_cuda2.h"

#include "../cuda/psc_cuda.h"

#define F3_CACHE F3_CUDA2

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_interpolate.c"
#include "../psc_push_particles/inc_push.c"
#include "../psc_push_particles/inc_curr.c"
#include "../psc_push_particles/inc_step.c"

