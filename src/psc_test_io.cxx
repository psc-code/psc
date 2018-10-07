
#include "dim.hxx"
#include "grid.hxx"
#include "psc_fields_single.h"

#include <mrc_io.hxx>

#include <string>

// ======================================================================
// PscTestIo

struct PscTestIo
{
  static void run()
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
    auto grid = Grid_t{grid_domain, grid_bc, kinds, norm, dt};

    mpi_printf(MPI_COMM_WORLD, "***** Testing output\n");

    auto io_pfd = MrcIo{"pfd", "."};
    auto io = io_pfd.io_;

    mrc_io_open(io_pfd.io_, "w", grid.timestep(), grid.timestep() * grid.dt);
    
    // save some basic info about the run in the output file
    struct mrc_obj *obj = mrc_obj_create(mrc_io_comm(io));
    mrc_obj_set_name(obj, "psc");
    mrc_obj_dict_add_int(obj, "timestep", grid.timestep());
    mrc_obj_dict_add_float(obj, "time", grid.timestep() * grid.dt);
    mrc_obj_dict_add_float(obj, "cc", grid.norm.cc);
    mrc_obj_dict_add_float(obj, "dt", grid.dt);
    mrc_obj_write(obj, io);
    mrc_obj_destroy(obj);
    
    mrc_fld* fld = grid.mrc_domain().m3_create();
    mrc_fld_set_name(fld, "e");
    mrc_fld_set_param_int(fld, "nr_ghosts", 0);
    mrc_fld_set_param_int(fld, "nr_comps", 2);
    mrc_fld_setup(fld);
    mrc_fld_set_comp_name(fld, 0, "ex");
    mrc_fld_set_comp_name(fld, 1, "ey");
    
    for (int p = 0; p < grid.n_patches(); p++) {
      mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
      mrc_fld_foreach(fld, i,j,k, 0,0) {
	MRC_M3(m3p, 0, i,j,k) = i;
	MRC_M3(m3p, 1, i,j,k) = j;
      } mrc_fld_foreach_end;
      mrc_fld_patch_put(fld);
    }
    
    mrc_fld_write(fld, io);
    mrc_fld_destroy(fld);

    mrc_io_close(io);

    mpi_printf(MPI_COMM_WORLD, "***** Testing output done\n");
  }
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  PscTestIo::run();

  libmrc_params_finalize();
  MPI_Finalize();
  return 0;
}
