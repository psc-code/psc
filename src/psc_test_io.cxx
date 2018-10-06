
#include "dim.hxx"
#include "grid.hxx"
#include "psc_fields_single.h"

#include <mrc_io.hxx>

#include <string>

// ======================================================================
// PscTestIo

struct PscTestIo
{
  // ----------------------------------------------------------------------
  // ctor
  
  PscTestIo()
  {
    mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

    // -- setup particle kinds
    Grid_t::Kinds kinds = {};
    
    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
#if 0
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 8, 16, 48 }; // division into patches
#else
    Int3 gdims = { 20, 20, 80}; // global number of grid points
    Int3 np = { 2, 2, 8 }; // division into patches
#endif
    
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 5;
    
    double dt = .99;
    auto norm = Grid_t::Normalization{norm_params};
    grid_ = new Grid_t{grid_domain, grid_bc, kinds, norm, dt};

    mpi_printf(MPI_COMM_WORLD, "***** Testing output\n");

    Int3 rn = {};
    Int3 rx = {1000000, 1000000, 100000};

    auto io_pfd = MrcIo{"pfd", "."};
    io_pfd.open(grid(), rn, rx);

    auto mres = MfieldsSingle{grid(), 2, {2, 2, 2}};
    mres.write_as_mrc_fld(io_pfd.io_, "e", {"ex", "ey"});

    io_pfd.close();

    mpi_printf(MPI_COMM_WORLD, "***** Testing output done\n");
  }

  const Grid_t& grid() { return *grid_; }

protected:
  Grid_t* grid_;
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  auto psc = PscTestIo{};
  
  MPI_Finalize();
  return 0;
}
