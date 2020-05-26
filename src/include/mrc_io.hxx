
#pragma once

#include <mrc_io.h>

#include <cstring>

// ======================================================================
// MrcIo

namespace MrcIo
{

  // static version so it can be used elsewhere without MrcIo wrapper
template <typename Mfields>
static void write_mflds(mrc_io* io, const Mfields& _mflds, const Grid_t& grid,
                        const std::string& name,
                        const std::vector<std::string>& comp_names)
{
  auto& mflds = const_cast<Mfields&>(_mflds);
  int n_comps = comp_names.size();
  // FIXME, should generally equal the # of component in mflds,
  // but this way allows us to write fewer components, useful to hack around
  // 16-bit vpic material ids, stored together as AOS with floats...

  mrc_fld* fld = grid.mrc_domain().m3_create();
  mrc_fld_set_name(fld, name.c_str());
  mrc_fld_set_param_int(fld, "nr_ghosts", 0);
  mrc_fld_set_param_int(fld, "nr_comps", n_comps);
  mrc_fld_setup(fld);
  for (int m = 0; m < n_comps; m++) {
    mrc_fld_set_comp_name(fld, m, comp_names[m].c_str());
  }

  for (int p = 0; p < mflds.n_patches(); p++) {
    mrc_fld_patch* m3p = mrc_fld_patch_get(fld, p);
    mrc_fld_foreach(fld, i, j, k, 0, 0)
    {
      for (int m = 0; m < n_comps; m++) {
        MRC_M3(m3p, m, i, j, k) = mflds[p](m, i, j, k);
      }
    }
    mrc_fld_foreach_end;
    mrc_fld_patch_put(fld);
  }

  mrc_fld_write(fld, io);
  mrc_fld_destroy(fld);
}

} // namespace MrcIo
