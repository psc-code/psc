
#include <psc_config.h>

#include <psc.h>

#include <mrc_profile.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>


#include "fields_item.hxx"

#include <mrc_io.hxx>

#include <balance.hxx>
#include <particles.hxx>
#include <fields3d.hxx>
#include <push_particles.hxx>
#include <push_fields.hxx>
#include <sort.hxx>
#include <collision.hxx>
#include <bnd_particles.hxx>
#include <bnd.hxx>
#include <bnd_fields.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>

#include "../libpsc/psc_checks/checks_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
#endif

#include "psc_config.hxx"

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// EDIT to change order / floating point type / cuda / 2d/3d
using dim_t = dim_xyz;
using PscConfig = PscConfig1vbecSingle<dim_t>;

// ======================================================================
// PscTestIo

struct PscTestIo
{
  using Mparticles = PscConfig::Mparticles_t;
  using MfieldsState = PscConfig::MfieldsState;
  using Mfields = MfieldsSingle;
  using DIM = PscConfig::dim_t;

  // ----------------------------------------------------------------------
  // ctor
  
  PscTestIo()
    : grid_{ggrid}
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    BB_ = 0.;
    
    // -- setup particle kinds
    // last population ("e") is neutralizing
    // FIXME, hardcoded mass ratio 100
    Grid_t::Kinds kinds = {{1., 100., "i"}, { -1., 1., "e"}};
    
    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 40, 80, 4 }; // division into patches

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

    mflds_.reset(new MfieldsState{grid()});
    mprts_.reset(new Mparticles{grid()});

    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
	switch (m) {
	case HY: return BB_;
	default: return 0.;
	}
      });
  }

  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");

    std::string name = "e";

    auto item_ = psc_output_fields_item_create(grid().comm());
    psc_output_fields_item_set_type(item_, name.c_str());
    psc_output_fields_item_setup(item_);

    PscFieldsItemBase{item_}(*mflds_, *mprts_);
    
    mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n");

    Int3 rn = {};
    Int3 rx = {1000000, 1000000, 100000};

    auto io_pfd = MrcIo{"pfd", "."};
    io_pfd.open(grid(), rn, rx);

    auto comp_names = PscFieldsItemBase{item_}->comp_names();
    auto& pfd = PscFieldsItemBase{item_}->mres();
    auto mres = Mfields{grid(), 3, grid().ibn};
    mres.write_as_mrc_fld(io_pfd.io_, name, comp_names);

    io_pfd.close();

    psc_output_fields_item_destroy(item_);

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  const Grid_t& grid() { return *grid_; }

private:
  double BB_;

protected:
  Grid_t*& grid_;
  std::unique_ptr<MfieldsState> mflds_;
  std::unique_ptr<Mparticles> mprts_;

  Int3 ibn = {2, 2, 2}; // FIXME!!! need to factor in invar dims (but not in vpic...)
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);
  
  auto psc = new PscTestIo;

  psc->initialize();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
