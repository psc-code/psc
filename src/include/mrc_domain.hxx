
#pragma once

#include "grid/BC.h"
#include "grid/domain.hxx"

#include <mrc_domain_private.h>

#include <vector>

// ======================================================================
// MrcDomain

struct MrcDomain
{
  MrcDomain() = default;

  template <typename R>
  MrcDomain(const psc::grid::Domain<R>& grid_domain,
            const psc::grid::BC& grid_bc, int nr_patches)
  {
    domain_ = mrc_domain_create(MPI_COMM_WORLD);
    // create a very simple domain decomposition
    int bc[3] = {};
    for (int d = 0; d < 3; d++) {
      if (grid_bc.fld_lo[d] == BND_FLD_PERIODIC && grid_domain.gdims[d] > 1) {
        bc[d] = BC_PERIODIC;
      }
    }

    mrc_domain_set_type(domain_, "multi");
    mrc_domain_set_param_int3(domain_, "m", grid_domain.gdims);
    mrc_domain_set_param_int(domain_, "bcx", bc[0]);
    mrc_domain_set_param_int(domain_, "bcy", bc[1]);
    mrc_domain_set_param_int(domain_, "bcz", bc[2]);
    mrc_domain_set_param_int(domain_, "nr_patches", nr_patches);
    mrc_domain_set_param_int3(domain_, "np", grid_domain.np);

    struct mrc_crds* crds = mrc_domain_get_crds(domain_);
    mrc_crds_set_type(crds, "uniform");
    mrc_crds_set_param_int(crds, "sw", 2);
    mrc_crds_set_param_double3(crds, "l", grid_domain.corner);
    mrc_crds_set_param_double3(crds, "h",
                               grid_domain.corner + grid_domain.length);

    mrc_domain_set_from_options(domain_);
    mrc_domain_setup(domain_);

    // make sure that np isn't overridden on the command line
    Int3 np;
    mrc_domain_get_param_int3(domain_, "np", np);
    assert(np == Int3::fromPointer(grid_domain.np));
  }

  MrcDomain(const MrcDomain&) = delete;
  MrcDomain& operator=(const MrcDomain&) = delete;

  MrcDomain(MrcDomain&& o)
  {
    domain_ = o.domain_;
    o.domain_ = nullptr;
  }

  MrcDomain& operator=(MrcDomain&& o)
  {
    if (this != &o) {
      domain_ = o.domain_;
      o.domain_ = nullptr;
    }
    return *this;
  }

  ~MrcDomain()
  {
    if (domain_) {
      mrc_domain_destroy(domain_);
    }
  }

  MPI_Comm comm() const { return mrc_domain_comm(domain_); }
  void view() const { mrc_domain_view(domain_); }

  int nGlobalPatches() const
  {
    int n_global_patches;
    mrc_domain_get_nr_global_patches(domain_, &n_global_patches);
    return n_global_patches;
  }
  int nPatches() const
  {
    int n_patches;
    mrc_domain_get_patches(domain_, &n_patches);
    return n_patches;
  }
  mrc_patch_info globalPatchInfo(int p) const
  {
    mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain_, p, &info);
    return info;
  }
  mrc_patch_info localPatchInfo(int p) const
  {
    mrc_patch_info info;
    mrc_domain_get_local_patch_info(domain_, p, &info);
    return info;
  }
  mrc_patch_info levelIdx3PatchInfo(int level, int idx[3]) const
  {
    mrc_patch_info info;
    mrc_domain_get_level_idx3_patch_info(domain_, level, idx, &info);
    return info;
  }
  const mrc_patch* getPatches(int* n_patches) const
  {
    return mrc_domain_get_patches(domain_, n_patches);
  }
  void neighborRankPatch(int p, int dir[3], int* nei_rank, int* nei_patch) const
  {
    mrc_domain_get_neighbor_rank_patch(domain_, p, dir, nei_rank, nei_patch);
  }

  std::vector<Int3> offs() const
  {
    std::vector<Int3> offs;

    int n_patches;
    auto patches = getPatches(&n_patches);
    assert(n_patches > 0);
    Int3 ldims = Int3::fromPointer(patches[0].ldims);
    offs.reserve(n_patches);
    for (int p = 0; p < n_patches; p++) {
      assert(ldims == Int3::fromPointer(patches[p].ldims));
      offs.push_back(Int3::fromPointer(patches[p].off));
    }

    return offs;
  }

  mrc_fld* m3_create() const { return mrc_domain_m3_create(domain_); }
  mrc_ddc* create_ddc() const { return mrc_domain_create_ddc(domain_); }

private:
  mrc_domain* domain_ = nullptr;
};
