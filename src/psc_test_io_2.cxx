
#include <psc.h>

#include <mrc_profile.h>

#include "fields_item.hxx"

#include <mrc_io.hxx>

#include "psc_fields_single.h"

// ======================================================================
// PscTestIo

struct PscTestIo
{
  using Mfields = MfieldsSingle;
  using dim = dim_xyz;

  // ----------------------------------------------------------------------
  // ctor
  
  PscTestIo()
    : grid_{ggrid}
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    // -- setup particle kinds
    // last population ("e") is neutralizing
    // FIXME, hardcoded mass ratio 100
    Grid_t::Kinds kinds = {{1., 100., "i"}, { -1., 1., "e"}};
    
    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
#if 0
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches
#else
    Int3 gdims = { 40, 10, 20}; // global number of grid points
    Int3 np = { 4, 1, 2 }; // division into patches
#endif

    if (dim::InvarX::value) { ibn[0] = 0; } // FIXME, wrong place, not for VPIC...
    
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 5;

    double dt = .99;
    auto coeff = Grid_t::Normalization{norm_params};
    grid_ = Grid_t::psc_make_grid(grid_domain, grid_bc, kinds, coeff, dt, ibn);
  }

  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    // initial output / stats
    mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n");

    Int3 rn = {};
    Int3 rx = {1000000, 1000000, 100000};

    auto io_pfd = MrcIo{"pfd", "."};
    io_pfd.open(grid(), rn, rx);

    auto mres = Mfields{grid(), 3, grid().ibn};
    mres.write_as_mrc_fld(io_pfd.io_, "e", {"ex", "ey", "ez"});

    io_pfd.close();

    MHERE;
    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  const Grid_t& grid() { return *grid_; }

protected:
  Grid_t*& grid_;
  Int3 ibn = {2, 2, 2}; // FIXME!!! need to factor in invar dims (but not in vpic...)
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);

  {
    auto psc = PscTestIo{};
    psc.initialize();
  }

  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
