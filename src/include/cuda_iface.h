
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_particle_buf_cuda.h"
#include "psc_fields_single.h"

#include "cuda_mparticles.h"

// ----------------------------------------------------------------------
// cuda_base

void cuda_base_init(void);

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

struct cuda_mfields *cuda_mfields_create(void);
void cuda_mfields_destroy(struct cuda_mfields *cmflds);
void cuda_mfields_ctor(struct cuda_mfields *cmflds, mrc_json_t json);
void cuda_mfields_dtor(struct cuda_mfields *cmflds);
mrc_json_t cuda_mfields_to_json(struct cuda_mfields *cmflds);
void cuda_mfields_dump(struct cuda_mfields *cmflds, const char *filename);

fields_single_t cuda_mfields_get_host_fields(struct cuda_mfields *cmflds);
void cuda_mfields_copy_to_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me);
void cuda_mfields_copy_from_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me);
void cuda_mfields_axpy_comp_yz(struct cuda_mfields *y, int ym, float a, struct cuda_mfields *x, int xm);
void cuda_mfields_zero_comp_yz(struct cuda_mfields *x, int xm);

void cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p);

void cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt);
void cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt);
void cuda_marder_correct_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3], int lz[3], int rz[3]);

// ----------------------------------------------------------------------
// cuda_push_mprts_yz

void cuda_push_mprts_yz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds,
			int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);

// ----------------------------------------------------------------------
// cuda_moments

void cuda_moments_yz_rho_1st_nc(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres);
void cuda_moments_yz_n_1st(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres);

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
void cuda_heating_run_foil(struct cuda_mparticles *cmprts);

// FIXME, mv elsewhere
#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
