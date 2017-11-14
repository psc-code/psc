
#ifndef BK_MPARTICLES_IFACE_H
#define BK_MPARTICLES_IFACE_H

#include "psc_particle_single_by_kind.h"

#ifdef __cplusplus
#include "bk_mparticles.h"
typedef std::vector<particle_single_by_kind_t> particle_single_by_kind_buf_t;
typedef mparticles<particle_single_by_kind_buf_t> bk_mparticles;
#else
typedef struct bk_mparticles bk_mparticles;
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

bk_mparticles *bk_mparticles_new(int n_patches);
void bk_mparticles_delete(bk_mparticles *bkmprts);
void bk_mparticles_reserve_all(bk_mparticles *bkmprts, int n_prts_by_patch[]);
void bk_mparticles_resize_all(bk_mparticles *bkmprts, int n_prts_by_patch[]);
void bk_mparticles_size_all(bk_mparticles *bkmprts, int n_prts_by_patch[]);
int bk_mparticles_n_prts(bk_mparticles *bkmprts);
particle_single_by_kind_t *bk_mparticles_at_ptr(bk_mparticles *bkmprts, int p, int n);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif


