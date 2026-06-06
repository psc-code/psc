
#include "io_common.h"
#include <kg/io.h>

#include "fields.hxx"

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

    writer.put("ibn", mflds_cuda.ibn(), launch);

    auto&& h_mflds = gt::host_mirror(mflds_cuda.storage());
    gt::copy(mflds_cuda.storage(), h_mflds);

    std::vector<Int3> patchOffsets;
    for (int p = 0; p < mflds_cuda.n_patches(); p++) {
      patchOffsets.push_back(grid.patches[p].off);
    }
    write_4d(writer, grid.ldims, grid.domain.gdims, patchOffsets, h_mflds);
  }

  void get(kg::io::Engine& reader, Mfields& mflds_cuda,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    const auto& grid = mflds_cuda.grid();
    auto n_comps = mflds_cuda.n_comps();
    auto n_patches = mflds_cuda.n_patches();
    auto&& mflds = gt::host_mirror(mflds_cuda.storage());
    auto h_mflds = gt::empty<typename Mfields::real_t>(gt::shape(
      grid.ldims[0], grid.ldims[1], grid.ldims[2], n_comps, n_patches));
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

    psc::mflds::interior(grid, mflds) = h_mflds;
    gt::copy(mflds, mflds_cuda.storage());
  }
};
