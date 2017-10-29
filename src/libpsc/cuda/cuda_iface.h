
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include "mrc_json.h"
#include "psc_particle_buf_cuda.h"
#include "psc_fields_single.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

// ----------------------------------------------------------------------
// cuda_base

void cuda_base_init(void);

// ----------------------------------------------------------------------
// double3 / float_4

typedef double double_3[3];
typedef float float_4[4];

// ----------------------------------------------------------------------
// cuda_mparticles_prt

struct cuda_mparticles_prt {
  float xi[3];
  float pxi[3];
  int kind;
  float qni_wni;
};

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles;

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_ctor(struct cuda_mparticles *cuda_mprts, mrc_json_t json);
void cuda_mparticles_reserve_all(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch);
void cuda_mparticles_dump(struct cuda_mparticles *cuda_mprts);
void cuda_mparticles_dump_by_patch(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch);
void cuda_mparticles_setup_internals(struct cuda_mparticles *cmprts);
unsigned int cuda_mparticles_get_n_prts(struct cuda_mparticles *cmprts);
void cuda_mparticles_get_size_all(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch);
void cuda_mparticles_resize_all(struct cuda_mparticles *cmprts, const unsigned int *n_prts_by_patch);
void cuda_mparticles_set_particles(struct cuda_mparticles *cmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx),
				   void *ctx);
void cuda_mparticles_get_particles(struct cuda_mparticles *cmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct cuda_mparticles_prt *, int, void *),
				   void *ctx);
void cuda_mparticles_to_device(struct cuda_mparticles *cmprts, float_4 *xi4, float_4 *pxi4,
			       unsigned int n_prts, unsigned int off);
void cuda_mparticles_from_device(struct cuda_mparticles *cmprts, float_4 *xi4, float_4 *pxi4,
				 unsigned int n_prts, unsigned int off);
  
void cuda_mparticles_inject(struct cuda_mparticles *cmprts, struct cuda_mparticles_prt *buf,
			    unsigned int *buf_n_by_patch);

const particle_cuda_real_t *cuda_mparticles_patch_get_b_dxi(struct cuda_mparticles *cmprts, int p);
const int *cuda_mparticles_patch_get_b_mx(struct cuda_mparticles *cmprts, int p);

psc_particle_cuda_buf_t *cuda_mparticles_bnd_get_buffer(struct cuda_mparticles *cmprts, int p);
void cuda_mparticles_bnd_prep(struct cuda_mparticles *cmprts);
void cuda_mparticles_bnd_post(struct cuda_mparticles *cmprts);

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

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
