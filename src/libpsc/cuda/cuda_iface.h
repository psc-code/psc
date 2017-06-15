
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

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
// cuda_mparticles

struct cuda_mparticles;

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				     const struct cuda_domain_info *info);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
