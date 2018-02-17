
#ifndef BK_MPARTICLES_IFACE_H
#define BK_MPARTICLES_IFACE_H

#ifdef __cplusplus

#include "bk_mparticles.h"

struct particle_single_by_kind {
  typedef float real_type;
  
  real_type dx[3];
  int i;
  real_type ux[3];
  real_type w;
  int kind;
};

typedef std::vector<particle_single_by_kind> particle_single_by_kind_buf;
typedef mparticles_<particle_single_by_kind_buf> bk_mparticles;

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
void bk_mparticles_reserve_all(bk_mparticles *bkmprts, const uint n_prts_by_patch[]);
void bk_mparticles_resize_all(bk_mparticles *bkmprts, const uint n_prts_by_patch[]);
void bk_mparticles_size_all(bk_mparticles *bkmprts, uint n_prts_by_patch[]);
int bk_mparticles_n_prts(bk_mparticles *bkmprts);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif


