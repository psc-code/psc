
#include <psc.h>

#include <mrc_profile.h>

#include "fields_item.hxx"

#include <mrc_io.hxx>

#include "psc_fields_single.h"

// ======================================================================
// PscTestIo

struct PscTestIo
{
  static void run()
  {
    auto comm = MPI_COMM_WORLD;

    mpi_printf(comm, "*** Setting up...\n");

    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
#if 0
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches
#else
    Int3 gdims = { 40, 10, 20}; // global number of grid points
    Int3 np = { 4, 1, 2 }; // division into patches
#endif

    struct mrc_domain *domain = mrc_domain_create(comm);
    
    mrc_domain_set_type(domain, "multi");
    mrc_domain_set_param_int3(domain, "m", gdims);
    mrc_domain_set_param_int3(domain, "np", np);
    
    struct mrc_crds *crds = mrc_domain_get_crds(domain);
    mrc_crds_set_type(crds, "uniform");
    mrc_crds_set_param_int(crds, "sw", 2);
    mrc_crds_set_param_double3(crds, "l", -.5 * LL);
    mrc_crds_set_param_double3(crds, "h", LL);
    
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    
    int n_patches;
    mrc_patch* patches = mrc_domain_get_patches(domain, &n_patches);

    mpi_printf(MPI_COMM_WORLD, "*** Writing output\n");

    mrc_io* io = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io, "basename", "pfd");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
    mrc_io_view(io);

    mrc_io_open(io, "w", 0, 0.);
    
    mrc_fld* fld = mrc_domain_m3_create(domain);
    mrc_fld_set_name(fld, "e");
    mrc_fld_set_param_int(fld, "nr_ghosts", 0);
    mrc_fld_set_param_int(fld, "nr_comps", 3);
    mrc_fld_setup(fld);
    mrc_fld_set_comp_name(fld, 0, "ex");
    mrc_fld_set_comp_name(fld, 1, "ey");
    mrc_fld_set_comp_name(fld, 2, "ez");
    
    for (int p = 0; p < n_patches; p++) {
      mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
      mrc_fld_foreach(fld, i,j,k, 0,0) {
	int ii = i + patches[p].off[0];
	int jj = j + patches[p].off[1];
	int kk = k + patches[p].off[2];
	for (int m = 0; m < mrc_fld_nr_comps(fld); m++) {
	  MRC_M3(m3p, m , i,j,k) = ii + 100 * jj + 10000 * kk;
	}
      } mrc_fld_foreach_end;
      mrc_fld_patch_put(fld);
    }
    
    mrc_fld_write(fld, io);
    mrc_fld_destroy(fld);

    mrc_io_close(io);
    mrc_io_destroy(io);

    mpi_printf(mrc_domain_comm(domain), "*** Writing output done.\n");
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
