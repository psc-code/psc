
#include "kg/io.h"

#include "psc_fields_c.h"

inline kg::io::Dims makeDims(int m, const Int3& dims)
{
  return kg::io::Dims{static_cast<size_t>(m), static_cast<size_t>(dims[2]),
                      static_cast<size_t>(dims[1]),
                      static_cast<size_t>(dims[0])};
}

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
    const Grid_t& grid = mflds.grid();

    writer.put("ib", mflds.box().ib(), launch);
    writer.put("im", mflds.box().im(), launch);

    auto n_comps = mflds.n_comps();
    auto shape = makeDims(n_comps, grid.domain.gdims);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto start = makeDims(0, grid.patches[p].off);
      auto count = makeDims(n_comps, grid.ldims);
      auto ib = makeDims(0, -mflds.box().ib());
      auto im = makeDims(n_comps, mflds.box().im());
      writer.putVariable(const_cast<Mfields&>(mflds)[p].data(), launch,
                         shape, {start, count}, {ib, im}); // FIXME cast
    }
  }

  void get(kg::io::Engine& reader, Mfields& mflds,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    const Grid_t& grid = mflds.grid();

    // FIXME, should just check for consistency? (# ghosts might differ, too)
    // reader.get("ib", mflds.ib, launch);
    // reader.get("im", mflds.im, launch);

    auto& gdims = grid.domain.gdims;
    auto n_comps = mflds.n_comps();
    auto shape = makeDims(n_comps, grid.domain.gdims);
    assert(reader.variableShape<DataType>() == shape);
    for (int p = 0; p < mflds.n_patches(); p++) {
      auto start = makeDims(0, grid.patches[p].off);
      auto count = makeDims(n_comps, grid.ldims);
      auto ib = makeDims(0, -mflds.box().ib());
      auto im = makeDims(n_comps, mflds.box().im());
      reader.getVariable(mflds[p].data(), launch, {start, count}, {ib, im});
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
