
#include <psc_config.h>

#include <psc.h>

#include <mrc_profile.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>
#include <output_fields_c.hxx>

// ======================================================================
// Psc

template<typename PscConfig>
struct Psc
{
  using Mparticles_t = typename PscConfig::Mparticles_t;
  using MfieldsState = typename PscConfig::MfieldsState;

  // ----------------------------------------------------------------------
  // ctor

  Psc()
    : grid_{ggrid}
  {}

  // ----------------------------------------------------------------------
  // define_grid

  void define_grid(Grid_t::Domain& domain, GridBc& bc, Grid_t::Kinds& kinds,
		   double dt, Grid_t::NormalizationParams& norm_params)
  {
    auto coeff = Grid_t::Normalization{norm_params};
    grid_ = Grid_t::psc_make_grid(domain, bc, kinds, coeff, dt, ibn);
  }
  
  // ----------------------------------------------------------------------
  // define_field_array

  void define_field_array(double damp = 0.)
  {
    mflds_.reset(new MfieldsState{grid()});
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");
    diagnostics();

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

private:

  // ----------------------------------------------------------------------
  // diagnostics

  virtual void diagnostics()
  {
    (*outf_)(*mflds_, *mprts_);
  }

public:
  const Grid_t& grid() { return *grid_; }

protected:
  Grid_t*& grid_;

  std::unique_ptr<MfieldsState> mflds_;
  std::unique_ptr<Mparticles_t> mprts_;

  std::unique_ptr<OutputFieldsC> outf_;

  Int3 ibn = {2, 2, 2}; // FIXME!!! need to factor in invar dims (but not in vpic...)
};

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
//using PscConfig = PscConfig1vbecCuda<dim_t>;

// ======================================================================
// PscTestIo

struct PscTestIo : Psc<PscConfig>
{
  using DIM = PscConfig::dim_t;

  // ----------------------------------------------------------------------
  // ctor
  
  PscTestIo()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    BB_ = 0.;
    Zi_ = 1.;
    
    // -- setup particle kinds
    // last population ("e") is neutralizing
    // FIXME, hardcoded mass ratio 100
    Grid_t::Kinds kinds = {{Zi_, 100.*Zi_, "i"}, { -1., 1., "e"}};
    
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
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- output fields
    OutputFieldsCParams outf_params;
    outf_params.output_fields = "e";
    outf_params.pfield_step = 100;
    outf_.reset(new OutputFieldsC{grid(), outf_params});

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

private:
  double BB_;
  double Zi_;
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
