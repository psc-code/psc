
#include "bk_mparticles_iface.h"

// ======================================================================
// bk_mparticles C wrappers

bk_mparticles *bk_mparticles_new(int n_patches)
{
  return new bk_mparticles(n_patches);
}

void bk_mparticles_delete(bk_mparticles *bkmprts)
{
  delete bkmprts;
}

void bk_mparticles_reserve_all(bk_mparticles *bkmprts, const uint n_prts_by_patch[])
{
  bkmprts->reserve_all(n_prts_by_patch);
}

void bk_mparticles_resize_all(bk_mparticles *bkmprts, const uint n_prts_by_patch[])
{
  bkmprts->resize_all(n_prts_by_patch);
}

void bk_mparticles_size_all(bk_mparticles *bkmprts, uint n_prts_by_patch[])
{
  bkmprts->size_all(n_prts_by_patch);
}

int bk_mparticles_n_prts(bk_mparticles *bkmprts)
{
  return bkmprts->n_prts();
}

