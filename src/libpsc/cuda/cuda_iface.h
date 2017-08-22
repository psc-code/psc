
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#include <stdbool.h>

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
// cuda_domain_info
  
struct cuda_domain_info {
  int n_patches;
  int ldims[3]; // number of cells per patch
  int bs[3];    // size of each block (a.k.a. super-cell)
  double dx[3]; // size of a single cell
  double_3 *xb_by_patch;
};
  
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
void cuda_mparticles_setup(struct cuda_mparticles *cmprts);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				     const struct cuda_domain_info *info);
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

void cuda_mparticles_bnd_prep(struct cuda_mparticles *cmprts);
void cuda_mparticles_bnd_post(struct cuda_mparticles *cmprts);

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

struct cuda_mfields *cuda_mfields_create(void);
void cuda_mfields_destroy(struct cuda_mfields *cmflds);

// ----------------------------------------------------------------------
// cuda_push_mprts_yz

void cuda_push_mprts_yz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds,
			int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);

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
