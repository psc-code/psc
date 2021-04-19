
#include "io_common.h"
#include <kg/io.h>

// ======================================================================
// Variable<MfieldsCuda>
//
// FIXME, consolidate with host mfields write

template <typename MF>
class kg::io::Descr<
  MF, typename std::enable_if<std::is_same<MF, MfieldsCuda>::value ||
                                std::is_same<MF, MfieldsStateCuda>::value,
                              void>::type>
{
public:
  using Mfields = MF;
  using DataType = typename Mfields::real_t;

  void put(kg::io::Engine& writer, const Mfields& mflds_cuda,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    const auto& grid = mflds_cuda.grid();
    auto mflds = hostMirror(mflds_cuda);
    copy(mflds_cuda, mflds);

    writer.put("ib", mflds.box().ib(), launch);
    writer.put("im", mflds.box().im(), launch);

    auto n_comps = mflds.n_comps();
    auto shape = makeDims(n_comps, grid.domain.gdims);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto start = makeDims(0, grid.patches[p].off);
      auto count = makeDims(n_comps, grid.ldims);
      auto ib = makeDims(0, -mflds.box().ib());
      auto im = makeDims(n_comps, mflds.box().im());
      writer.putVariable(mflds[p].storage().data(), launch, shape,
                         {start, count}, {ib, im}); // FIXME cast
    }

    // host mirror will go away as this function returns, so need
    // to write
    writer.performPuts();
  }

  void get(kg::io::Engine& reader, Mfields& mflds_cuda,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    const auto& grid = mflds_cuda.grid();
    auto n_comps = mflds_cuda.n_comps();
    auto n_patches = mflds_cuda.n_patches();
    auto mflds = gt::host_mirror(mflds_cuda.gt());
    auto h_mflds = gt::empty(
      {grid.ldims[0], grid.ldims[1], grid.ldims[2], n_comps, n_patches});
    // FIXME, should just check for consistency? (# ghosts might differ, too)
    // reader.get("ib", mflds.ib, launch);
    // reader.get("im", mflds.im, launch);

    auto shape = makeDims(n_comps, grid.domain.gdims);
    assert(reader.variableShape<DataType>() == shape);
    for (int p = 0; p < n_patches; p++) {
      auto start = makeDims(0, grid.patches[p].off);
      auto count = makeDims(n_comps, grid.ldims);
      reader.getVariable(&h_mflds(0, 0, 0, 0, p), launch, {start, count});
    }
    reader.performGets();

    view_interior(mflds.gt(), mflds.ibn()) = h_mflds;
    gt::copy(mflds, mflds_cuda.gt());
  }
};
