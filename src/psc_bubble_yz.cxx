
#include <psc.h>
#include <psc.hxx>
#include "psc_config.hxx"
#include <psc_particles_single.h>
#include <psc_fields_single.h>
#ifdef USE_VPIC
#include "../libpsc/vpic/vpic_iface.h" // FIXME
#endif

#include "push_particles.hxx"
#include "push_fields.hxx"
#include "sort.hxx"
#include "collision.hxx"
#include "bnd_particles.hxx"
#include "balance.hxx"
#include "checks.hxx"
#include "marder.hxx"

#include "setup_particles.hxx"
#include "setup_fields.hxx"

#include <mrc_params.h>

#include <math.h>

// ======================================================================
// PscBubbleParams

struct PscBubbleParams
{
  double BB;
  double nnb;
  double nn0;
  double MMach;
  double LLn;
  double LLB;
  double LLz;
  double LLy;
  double TTe;
  double TTi;
  double MMi;
};

using PscConfig = PscConfig1vbecSingle<dim_yz>;

// ======================================================================
// PscBubble

struct PscBubble : Psc<PscConfig>, PscBubbleParams
{
  using DIM = PscConfig::dim_t;

  // ----------------------------------------------------------------------
  // ctor

  PscBubble()
  {
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    PscBubbleParams params;
    
    params.BB = .07;
    params.nnb = .1;
    params.nn0 = 1.;
    params.MMach = 3.;
    params.LLn = 200.;
    params.LLB = 200./6.;
    params.TTe = .02;
    params.TTi = .02;
    params.MMi = 100.;
    
    params.LLy = 2. * params.LLn;
    params.LLz = 3. * params.LLn;

    static_cast<PscBubbleParams&>(*this) = params;

    p_.nmax = 1000; //32000;
    p_.stats_every = 100;
  
    auto grid_domain = Grid_t::Domain{{1, 128, 512},
				      {params.LLn, params.LLy, params.LLz},
				      {0., -.5 * params.LLy, -.5 * params.LLz},
				      {1, 1, 4}};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto kinds = Grid_t::Kinds{{ -1., 1., "e"}, { 1., 100., "i" }};
    
    mpi_printf(comm, "lambda_D = %g\n", sqrt(params.TTe));
    
    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 100;

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles_t{grid()});

    // -- Balance
    balance_interval = 0;
    balance_.reset(new Balance_t{balance_interval, .1,});

    // -- Sort
    // FIXME, needs a way to make it gets set?
    sort_interval = 10;

    // -- Collision
    int collision_interval = 10;
    double collision_nu = .1;
    collision_.reset(new Collision_t{grid(), collision_interval, collision_nu});

    // -- Checks
    ChecksParams checks_params;
    checks_.reset(new Checks_t{grid(), comm, checks_params});

    // -- Marder correction
    double marder_diffusion = 0.9;
    int marder_loop = 3;
    bool marder_dump = false;
    marder_interval = 0*5;
    marder_.reset(new Marder_t(grid(), marder_diffusion, marder_loop, marder_dump));

    // -- output fields
    OutputFieldsCParams outf_params;
    outf_params.output_fields = "e,h,j,n_1st_single,v_1st_single";
    outf_params.pfield_step = 10;
    outf_.reset(new OutputFieldsC{grid(), outf_params});

    // --- partition particles and initial balancing
    mpi_printf(comm, "**** Partitioning...\n");
    auto n_prts_by_patch_old = setup_initial_partition();
    auto n_prts_by_patch_new = balance_->initial(n_prts_by_patch_old);
    // balance::initial does not rebalance particles, because the old way of doing this
    // does't even have the particle data structure created yet -- FIXME?
    mprts_->reset(grid());
    
    mpi_printf(comm, "**** Setting up particles...\n");
    setup_initial_particles(*mprts_, n_prts_by_patch_new);
    
    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);

    init();
  }

  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {
    double V0 = MMach * sqrt(TTe / MMi);
    
    double r1 = sqrt(sqr(crd[2]) + sqr(crd[1] + .5 * LLy));
    double r2 = sqrt(sqr(crd[2]) + sqr(crd[1] - .5 * LLy));
    
    npt.n = nnb;
    if (r1 < LLn) {
      npt.n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r1 / LLn));
      if (r1 > 0.0) {
	npt.p[2] += V0 * sin(M_PI * r1 / LLn) * crd[2] / r1;
	npt.p[1] += V0 * sin(M_PI * r1 / LLn) * (crd[1] + .5 * LLy) / r1;
      }
    }
    if (r2 < LLn) {
      npt.n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r2 / LLn));
      if (r2 > 0.0) {
	npt.p[2] += V0 * sin(M_PI * r2 / LLn) * crd[2] / r2;
	npt.p[1] += V0 * sin(M_PI * r2 / LLn) * (crd[1] - .5 * LLy) / r2;
      }
    }
    
    switch (kind) {
    case 0: // electrons
      // electron drift consistent with initial current
      if ((r1 <= LLn) && (r1 >= LLn - 2.*LLB)) {
	npt.p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r1)/(2.*LLB)) / npt.n;
      }
      if ((r2 <= LLn) && (r2 >= LLn - 2.*LLB)) {
	npt.p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r2)/(2.*LLB)) / npt.n;
      }
      
      npt.T[0] = TTe;
      npt.T[1] = TTe;
      npt.T[2] = TTe;
      break;
    case 1: // ions
      npt.T[0] = TTi;
      npt.T[1] = TTi;
      npt.T[2] = TTi;
      break;
    default:
      assert(0);
    }
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_partition
  
  std::vector<uint> setup_initial_partition()
  {
    SetupParticles<Mparticles_t> setup_particles;    
    return setup_particles.setup_partition(grid(), [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_particles
  
  void setup_initial_particles(Mparticles_t& mprts, std::vector<uint>& n_prts_by_patch)
  {
#if 0
    n_prts_by_patch[0] = 2;
    mprts.reserve_all(n_prts_by_patch.data());
    mprts.resize_all(n_prts_by_patch.data());

    for (int p = 0; p < mprts.n_patches(); p++) {
      mprintf("npp %d %d\n", p, n_prts_by_patch[p]);
      for (int n = 0; n < n_prts_by_patch[p]; n++) {
	auto &prt = mprts[p][n];
	prt.pxi = n;
	prt.kind_ = n % 2;
	prt.qni_wni_ = mprts.grid().kinds[prt.kind_].q;
      }
    };
#else
    SetupParticles<Mparticles_t> setup_particles;
    setup_particles.setup_particles(mprts, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
#endif
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    SetupFields<MfieldsState>::set(mflds, [&](int m, double crd[3]) {
	double z1 = crd[2];
	double y1 = crd[1] + .5 * LLy;
	double r1 = sqrt(sqr(z1) + sqr(y1));
	double z2 = crd[2];
	double y2 = crd[1] - .5 * LLy;
	double r2 = sqrt(sqr(z2) + sqr(y2));

	double rv = 0.;
	switch (m) {
	case HZ:
	  if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	    rv += - BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * y1 / r1;
	  }
	  if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	    rv += - BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * y2 / r2;
	  }
	  return rv;
	  
	case HY:
	  if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	    rv += BB * sin(M_PI * (LLn - r1)/(2.*LLB)) * z1 / r1;
	  }
	  if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	    rv += BB * sin(M_PI * (LLn - r2)/(2.*LLB)) * z2 / r2;
	  }
	  return rv;
	  
	case EX:
	  if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	    rv += MMach * sqrt(TTe/MMi) * BB *
	      sin(M_PI * (LLn - r1)/(2.*LLB)) * sin(M_PI * r1 / LLn);
	  }
	  if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	    rv += MMach * sqrt(TTe/MMi) * BB *
	      sin(M_PI * (LLn - r2)/(2.*LLB)) * sin(M_PI * r2 / LLn);
	  }
	  return rv;

	  // FIXME, JXI isn't really needed anymore (?)
	case JXI:
	  if ( (r1 < LLn) && (r1 > LLn - 2*LLB) ) {
	    rv += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r1)/(2.*LLB));
	  }
	  if ( (r2 < LLn) && (r2 > LLn - 2*LLB) ) {
	    rv += BB * M_PI/(2.*LLB) * cos(M_PI * (LLn - r2)/(2.*LLB));
	  }
	  return rv;
	  
	default:
	  return 0.;
	}
      });
  }
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);

  auto psc = new PscBubble;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
