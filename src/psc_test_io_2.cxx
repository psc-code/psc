
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

  static void run()
  {
    auto comm = MPI_COMM_WORLD;

    mpi_printf(comm, "*** Setting up...\n");

    // -- setup particle kinds
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
    auto& grid = *Grid_t::psc_make_grid(grid_domain, grid_bc, kinds, coeff, dt, {2, 2, 2});

    mpi_printf(MPI_COMM_WORLD, "*** Writing output\n");

    Int3 rn = {};
    Int3 rx = {1000000, 1000000, 100000};

    auto io_pfd = MrcIo{"pfd", "."};
    io_pfd.open(grid, rn, rx);

    auto mres = Mfields{grid, 3, grid.ibn};

    for (int p = 0; p < mres.n_patches(); p++) {
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	  int ii = i + grid.patches[p].off[0];
	  int jj = j + grid.patches[p].off[1];
	  int kk = k + grid.patches[p].off[2];
	  mres[p](0, i,j,k) = ii + 100 * jj + 10000 * kk;
	});
    }
    MrcIo::write_mflds(io_pfd.io_, mres, "e", {"ex", "ey", "ez"});

    io_pfd.close();

    mpi_printf(grid.comm(), "*** Writing output done.\n");
  }
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);

  PscTestIo::run();

  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
