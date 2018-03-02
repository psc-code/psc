
#include "psc_debug.h"

#define MAX_NR_KINDS (10)

#define PARTICLE_LOAD(prt, mprts_arr, n)	\
  prt = &mprts_arr[n]

#define PARTICLE_STORE(prt, mprts_arr, n) do {} while (0)

// ======================================================================

template<typename C>
struct PushOne
{
};
