
#ifndef CUDA_IFACE_H
#define CUDA_IFACE_H

#ifdef __cplusplus
extern "C" {
#endif

struct cuda_mparticles;

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);

#ifdef __cplusplus
}
#endif

#endif
