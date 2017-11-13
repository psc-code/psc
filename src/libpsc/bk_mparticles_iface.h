
#ifndef BK_MPARTICLES_IFACE_H
#define BK_MPARTICLES_IFACE_H

#include "psc_particle_single_by_kind.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

struct bk_mparticles;

struct bk_mparticles *bk_mparticles_create();
void bk_mparticles_destroy(struct bk_mparticles *bkmprts);
void bk_mparticles_ctor(struct bk_mparticles *bkmprts, int n_patches);
void bk_mparticles_dtor(struct bk_mparticles *bkmprts);
void bk_mparticles_reserve_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[]);
void bk_mparticles_resize_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[]);
void bk_mparticles_get_size_all(struct bk_mparticles *bkmprts, int n_prts_by_patch[]);
int bk_mparticles_get_n_prts(struct bk_mparticles *bkmprts);
particle_single_by_kind_t *bk_mparticles_at_ptr(struct bk_mparticles *bkmprts, int p, int n);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif


