
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
// cuda_mfields

struct cuda_mfields;

void cuda_push_fields_E_xyz(struct cuda_mfields* cmflds, float dt);
void cuda_push_fields_H_xyz(struct cuda_mfields* cmflds, float dt);

// ----------------------------------------------------------------------
// cuda_moments

template <typename BS>
struct cuda_mparticles;

template <typename BS>
void cuda_moments_yz_rho_1st_nc(cuda_mparticles<BS>* cmprts,
                                struct cuda_mfields* cmres);

template <typename BS>
void cuda_moments_yz_n_1st(cuda_mparticles<BS>* cmprts,
                           struct cuda_mfields* cmres);

// FIXME, mv elsewhere
#define HERE                                                                   \
  printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
