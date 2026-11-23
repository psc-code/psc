#pragma once

#include "turbulence_base.hxx"

#include "io_common.h"
#include "kg/io.h"

#include <string>
#include <vector>

namespace psc
{
namespace turbulence
{

namespace detail
{

// A contiguous segment of a patch's index range along one dim, mapped into the
// periodic global domain.
struct PeriodicSegment
{
  int local_start;
  int global_start;
  int count;
};

// Splits [ib, ib + im) (patch-local, starting at global offset `off`) into
// segments that don't cross the global domain's periodic boundary.
inline std::vector<PeriodicSegment> periodic_segments(int ib, int im, int off,
                                                      int gdim)
{
  std::vector<PeriodicSegment> segments;
  for (int i = ib; i < ib + im;) {
    int global_start = ((off + i) % gdim + gdim) % gdim;
    int count = std::min(ib + im - i, gdim - global_start);
    segments.push_back({i, global_start, count});
    i += count;
  }
  return segments;
}

// For kg::io::Engine::put<Var> / get<Var>. Writes H as one global array
// without ghosts, so the file doesn't depend on the decomposition. Reads it
// back with ghosts, which hold the periodic image of the domain.
template <typename MfieldsState>
struct PeriodicHVar
{
  using real_t = typename MfieldsState::real_t;
  static const int N_COMPS = 3;

  void put(kg::io::Engine& writer, const MfieldsState& mflds_const)
  {
    auto& mflds = const_cast<MfieldsState&>(mflds_const);
    const auto& grid = mflds.grid();
    Int3 ldims = grid.ldims;

    for (int p = 0; p < mflds.n_patches(); p++) {
      auto f = make_Fields3d<dim_xyz>(mflds[p]);
      std::vector<real_t> buf;
      buf.reserve(N_COMPS * ldims.prod());
      for (int m = 0; m < N_COMPS; m++) {
        grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
          buf.push_back(f(HX + m, jx, jy, jz));
        });
      }
      writer.putVariable(
        buf.data(), kg::io::Mode::Blocking,
        makeDims(N_COMPS, grid.domain.gdims),
        {makeDims(0, grid.patches[p].off), makeDims(N_COMPS, ldims)}, {});
    }
    writer.performPuts();
  }

  void get(kg::io::Engine& reader, MfieldsState& mflds)
  {
    const auto& grid = mflds.grid();
    Int3 gdims = grid.domain.gdims;
    Int3 ibn = mflds.ibn();

    if (reader.variableShape<real_t>() != makeDims(N_COMPS, gdims)) {
      LOG_ERROR("turbulence file doesn't match the grid's gdims\n");
    }

    for (int p = 0; p < mflds.n_patches(); p++) {
      auto f = make_Fields3d<dim_xyz>(mflds[p]);
      Int3 off = grid.patches[p].off;
      std::vector<PeriodicSegment> segments[3];
      for (int d = 0; d < 3; d++) {
        segments[d] = periodic_segments(-ibn[d], grid.ldims[d] + 2 * ibn[d],
                                        off[d], gdims[d]);
      }

      for (auto& sz : segments[2]) {
        for (auto& sy : segments[1]) {
          for (auto& sx : segments[0]) {
            Int3 local_start{sx.local_start, sy.local_start, sz.local_start};
            Int3 global_start{sx.global_start, sy.global_start,
                              sz.global_start};
            Int3 count{sx.count, sy.count, sz.count};

            std::vector<real_t> buf(N_COMPS * count.prod());
            reader.getVariable(
              buf.data(), kg::io::Mode::Blocking,
              {makeDims(0, global_start), makeDims(N_COMPS, count)}, {});

            int i = 0;
            for (int m = 0; m < N_COMPS; m++) {
              for (int jz = 0; jz < count[2]; jz++) {
                for (int jy = 0; jy < count[1]; jy++) {
                  for (int jx = 0; jx < count[0]; jx++) {
                    f(HX + m, local_start[0] + jx, local_start[1] + jy,
                      local_start[2] + jz) = buf[i++];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

} // namespace detail

#ifdef PSC_HAVE_ADIOS2

/// @brief Writes the turbulence in `mflds`' H to `filename` (an ADIOS2 .bp).
///
/// Only the interior is written; the file can be read back on any
/// decomposition of the same global grid.
template <typename MfieldsState>
void write_turbulence(const std::string& filename, MfieldsState& mflds)
{
  auto io = kg::io::IOAdios2{};
  auto writer =
    io.open(filename, kg::io::Mode::Write, mflds.grid().comm(), "turbulence");
  writer.template put<detail::PeriodicHVar>("h", mflds);
  writer.close();
}

/// @brief Reads turbulence written by write_turbulence into `mflds`' H.
///
/// Ghosts are filled as if the domain were periodic, which the generated
/// turbulence is, regardless of the grid's boundary conditions. Other
/// components are left alone.
template <typename MfieldsState>
void read_turbulence(const std::string& filename, MfieldsState& mflds)
{
  auto io = kg::io::IOAdios2{};
  auto reader =
    io.open(filename, kg::io::Mode::Read, mflds.grid().comm(), "turbulence");
  reader.beginStep(kg::io::StepMode::Read);
  reader.template get<detail::PeriodicHVar>("h", mflds);
  reader.endStep();
  reader.close();
}

#endif

} // namespace turbulence
} // namespace psc
