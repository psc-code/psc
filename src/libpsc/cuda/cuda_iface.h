
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_fields_single.h"

#include "psc_particles_cuda.h"
#include <particles_simple.hxx> // FIXME?
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
// cuda_base

void cuda_base_init(void);

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

void cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p);

void cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt);
void cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt);
void cuda_marder_correct_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3], int lz[3], int rz[3]);

// ----------------------------------------------------------------------
// cuda_push_mprts

template<typename BS>
struct cuda_mparticles;

void cuda_push_mprts_yz(cuda_mparticles<BS144>* cmprts, struct cuda_mfields *cmflds,
			const int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);

void cuda_push_mprts_xyz(cuda_mparticles<BS144>*cmprts, struct cuda_mfields *cmflds);

// ----------------------------------------------------------------------
// cuda_moments

void cuda_moments_yz_rho_1st_nc(cuda_mparticles<BS144>* cmprts, struct cuda_mfields *cmres);
void cuda_moments_yz_n_1st(cuda_mparticles<BS144>* cmprts, struct cuda_mfields *cmres);

// ----------------------------------------------------------------------
// cuda_heating_run_foil

struct cuda_heating_foil {
  // params
  float zl;
  float zh;
  float xc;
  float yc;
  float rH;
  float T;
  float Mi;
  int kind;

  // state (FIXME, shouldn't be part of the interface)
  float fac;
  float heating_dt;
};

void cuda_heating_setup_foil(struct cuda_heating_foil *foil);
void cuda_heating_run_foil(cuda_mparticles<BS144>* cmprts);

// FIXME, mv elsewhere
#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
