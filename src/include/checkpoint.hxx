
#pragma once

#include "fields3d.inl"
#include "grid.hxx"
#include "grid.inl"
#include "particles_simple.inl"
#include <kg/io.h>

// ----------------------------------------------------------------------
// write_checkpoint
//

template <typename Mparticles, typename MfieldsState>
void write_checkpoint(const Grid_t& grid, Mparticles& mprts,
                      MfieldsState& mflds)
{
  mpi_printf(grid.comm(), "**** Writing checkpoint...\n");
#if defined(PSC_HAVE_ADIOS2) && !defined(VPIC)
  MPI_Barrier(grid.comm()); // not really necessary

  std::string filename =
    "checkpoint_" + std::to_string(grid.timestep()) + ".bp";

  auto io = kg::io::IOAdios2{};
  auto writer = io.open(filename, kg::io::Mode::Write);
  writer.put("grid", grid);
  writer.put("mprts", mprts);
  writer.put("mflds", mflds);
  writer.close();
#else
  std::cerr << "write_checkpoint not available without adios2" << std::endl;
  std::abort();
#endif
}

// ----------------------------------------------------------------------
// read_checkpoint
//

template <typename Mparticles, typename MfieldsState>
inline void read_checkpoint(const std::string& filename, Grid_t& grid,
                            Mparticles& mprts, MfieldsState& mflds)
{
  mpi_printf(grid.comm(), "**** Reading checkpoint...\n");
  MPI_Barrier(grid.comm()); // not really necessary

#ifdef PSC_HAVE_ADIOS2
  auto io = kg::io::IOAdios2{};
  auto reader = io.open(filename, kg::io::Mode::Read);
  reader.get("grid", grid);
  reader.get("mprts", mprts);
  reader.get("mflds", mflds);
  reader.close();

  // FIXME, when we read back a rebalanced grid, other existing fields will
  // still have their own parallel distribution, ie, things will go wrong
#else
  std::cerr << "write_checkpoint not available without adios2" << std::endl;
  std::abort();
#endif
}
