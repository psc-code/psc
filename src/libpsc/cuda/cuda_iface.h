
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
// cuda_domain_info
  
struct cuda_domain_info {
  int n_patches;
  int ldims[3]; // number of cells per patch
  int bs[3];    // size of each block (a.k.a. super-cell)
  double dx[3]; // size of a single cell
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
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				     const struct cuda_domain_info *info);
void cuda_mparticles_alloc(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch);
void cuda_mparticles_dealloc(struct cuda_mparticles *cmprts);
void cuda_mparticles_dump(struct cuda_mparticles *cuda_mprts);
void cuda_mparticles_sort_initial(struct cuda_mparticles *cmprts,
				  unsigned int *n_prts_by_patch);
void cuda_mparticles_set_particles(struct cuda_mparticles *cmprts, unsigned int n_prts, unsigned int off,
				   void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx),
				   void *ctx);
void cuda_mparticles_get_particles(struct cuda_mparticles *cmprts, unsigned int n_prts, unsigned int off,
				   void (*put_particle)(struct cuda_mparticles_prt *, int, void *),
				   void *ctx);

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields;

struct cuda_mfields *cuda_mfields_create(void);
void cuda_mfields_destroy(struct cuda_mfields *cmflds);

// ----------------------------------------------------------------------
// cuda_push_mprts_yz

void cuda_push_mprts_yz(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds,
			int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);


#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
