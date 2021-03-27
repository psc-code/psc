
#pragma once

#include "io_common.h"
#include "kg/io.h"

// ======================================================================
// Variable<Mfields>

template <typename R>
class kg::io::Descr<Mfields<R>>
{
public:
  using Mfields = ::Mfields<R>;
  using DataType = typename Mfields::real_t;

  void put(kg::io::Engine& writer, const Mfields& mflds,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("ib", mflds.box().ib(), launch);
    writer.put("im", mflds.box().im(), launch);

    auto n_comps = mflds.n_comps();
    auto shape = makeDims(n_comps, mflds.gdims());
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto start = makeDims(0, mflds.patchOffset(p));
      auto count = makeDims(n_comps, mflds.ldims());
      auto ib = makeDims(0, -mflds.box().ib());
      auto im = makeDims(n_comps, mflds.box().im());
      writer.putVariable(const_cast<Mfields&>(mflds)[p].storage().data(),
                         launch, shape, {start, count}, {ib, im}); // FIXME cast
    }
  }

  void get(kg::io::Engine& reader, Mfields& mflds,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    // FIXME, should just check for consistency? (# ghosts might differ, too)
    // reader.get("ib", mflds.ib, launch);
    // reader.get("im", mflds.im, launch);

    auto n_comps = mflds.n_comps();
    auto shape = makeDims(n_comps, mflds.gdims());
    assert(reader.variableShape<DataType>() == shape);
    // FIXME, working around adios2 bug with selection
    MfieldsSingle h_mflds(mflds.grid(), n_comps, {});
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto start = makeDims(0, mflds.patchOffset(p));
      auto count = makeDims(n_comps, mflds.ldims());
      // auto ib = makeDims(0, -mflds.box().ib());
      // auto im = makeDims(n_comps, mflds.box().im());
      reader.getVariable(h_mflds[p].storage().data(), launch, {start, count},
                         {}); //{ib, im});
    }
    reader.performGets();

    for (int p = 0; p < mflds.n_patches(); p++) {
      auto h_flds = make_Fields3d<dim_xyz>(h_mflds[p]);
      auto flds = make_Fields3d<dim_xyz>(mflds[p]);
      for (int m = 0; m < n_comps; m++) {
        h_mflds.Foreach_3d(0, 0, [&](int i, int j, int k) {
          flds(m, i, j, k) = h_flds(m, i, j, k);
        });
      }
    }
  }
};

template <typename MFields>
class kg::io::Descr<MfieldsStateFromMfields<MFields>>
{
public:
  using MfieldsState = MfieldsStateFromMfields<MFields>;
  using DataType = typename MfieldsState::real_t;

  void put(kg::io::Engine& writer, const MfieldsState& mflds,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("mflds", mflds.mflds_);
  }

  void get(kg::io::Engine& reader, MfieldsState& mflds,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    reader.get("mflds", mflds.mflds_);
  }
};
