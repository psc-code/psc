
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_fields_single.h"

#include <grid.hxx>

#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...)                                                           \
  do {                                                                         \
  } while (0)
#endif

// ----------------------------------------------------------------------
// float_3 etc

typedef float float_3[3];
typedef double double_3[3];

// ----------------------------------------------------------------------
// cuda_moments

template <typename BS>
struct cuda_mparticles;

// FIXME, mv elsewhere
#define HERE                                                                   \
  printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
