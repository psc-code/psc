
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_fields_single.h"

#include <grid.hxx>

#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

// ----------------------------------------------------------------------
// float_3 etc

typedef float float_3[3];
typedef double double_3[3];

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

void cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p);

void cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt);
void cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt);
void cuda_marder_correct_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3], int lz[3], int rz[3]);

void cuda_push_fields_E_xyz(struct cuda_mfields *cmflds, float dt);
void cuda_push_fields_H_xyz(struct cuda_mfields *cmflds, float dt);
void cuda_marder_correct_xyz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			     int p, float fac[3],
			     int lx[3], int rx[3], int ly[3], int ry[3], int lz[3], int rz[3]);

// ----------------------------------------------------------------------
// cuda_moments

template<typename BS>
struct cuda_mparticles;

template<typename BS>
void cuda_moments_yz_rho_1st_nc(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres);

template<typename BS>
void cuda_moments_yz_n_1st(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres);

// FIXME, mv elsewhere
#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
