
#include <psc.h>
#include <psc.hxx>
#include "psc_config.hxx"
#include <psc_particles_single.h>
#include <psc_fields_single.h>

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

#include "../src/libpsc/psc_output_fields/fields_item_moments_1st.hxx"

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

// ======================================================================
// Global parameters

namespace
{

// Parameters specific to this run. They don't really need to be collected in a struct,
// but maybe it's nice that they are
PscBubbleParams g;

// This is a set of generic PSC params (see include/psc.hxx),
// like number of steps to run, etc, which also should be set by the case
PscParams psc_params;

}

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using PscConfig = PscConfig1vbecSingle<dim_yz>;

// ======================================================================
// PscBubble

struct PscBubble : Psc<PscConfig>
{
  using DIM = PscConfig::dim_t;

  // ----------------------------------------------------------------------
  // ctor

  PscBubble(const PscParams& psc_params)
  {
    p_ = psc_params;
    
    auto comm = grid().comm();

    mpi_printf(comm, "*** Setting up...\n");

    auto grid_domain = Grid_t::Domain{{1, 128, 512},
				      {g.LLn, g.LLy, g.LLz},
				      {0., -.5 * g.LLy, -.5 * g.LLz},
				      {1, 1, 4}};
    
    auto grid_bc = psc::grid::BC{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto kinds = Grid_t::Kinds{{ -1., 1., "e"}, { 1., 100., "i" }};
    
    mpi_printf(comm, "lambda_D = %g\n", sqrt(g.TTe));
    
    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 100;

    double dt = p_.cfl * courant_length(grid_domain);
    define_grid(grid_domain, grid_bc, kinds, dt, norm_params);

    define_field_array();

    mprts_.reset(new Mparticles{grid()});

    // -- Balance
    p_.balance_interval = 0;
    balance_.reset(new Balance_t{p_.balance_interval, .1,});

    // -- Sort
    // FIXME, needs a way to make it gets set?
    p_.sort_interval = 10;

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
    p_.marder_interval = 0*5;
    marder_.reset(new Marder_t(grid(), marder_diffusion, marder_loop, marder_dump));

    // -- output fields
    OutputFieldsCParams outf_params;
    outf_params.pfield_step = 10;
    std::vector<std::unique_ptr<FieldsItemBase>> outf_items;
    outf_items.emplace_back(new FieldsItemFields<ItemLoopPatches<Item_e_cc>>(grid()));
    outf_items.emplace_back(new FieldsItemFields<ItemLoopPatches<Item_h_cc>>(grid()));
    outf_items.emplace_back(new FieldsItemFields<ItemLoopPatches<Item_j_cc>>(grid()));
    outf_items.emplace_back(new FieldsItemMoment<ItemMomentAddBnd<Moment_n_1st<Mparticles, MfieldsC>>>(grid()));
    outf_items.emplace_back(new FieldsItemMoment<ItemMomentAddBnd<Moment_v_1st<Mparticles, MfieldsC>>>(grid()));
    outf_items.emplace_back(new FieldsItemMoment<ItemMomentAddBnd<Moment_T_1st<Mparticles, MfieldsC>>>(grid()));
    outf_.reset(new OutputFieldsC{grid(), outf_params, std::move(outf_items)});

    // --- partition particles and initial balancing
    mpi_printf(comm, "**** Partitioning...\n");
    auto n_prts_by_patch = setup_initial_partition();
    balance_->initial(grid_, n_prts_by_patch);
    // balance::initial does not rebalance particles, because the old way of doing this
    // does't even have the particle data structure created yet -- FIXME?
    mprts_->reset(grid());
    
    mpi_printf(comm, "**** Setting up particles...\n");
    setup_initial_particles(*mprts_, n_prts_by_patch);
    
    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(*mflds_);

    init();
  }

  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {
    double V0 = g.MMach * sqrt(g.TTe / g.MMi);
    
    double r1 = sqrt(sqr(crd[2]) + sqr(crd[1] + .5 * g.LLy));
    double r2 = sqrt(sqr(crd[2]) + sqr(crd[1] - .5 * g.LLy));
    
    npt.n = g.nnb;
    if (r1 < g.LLn) {
      npt.n += (g.nn0 - g.nnb) * sqr(cos(M_PI / 2. * r1 / g.LLn));
      if (r1 > 0.0) {
	npt.p[2] += V0 * sin(M_PI * r1 / g.LLn) * crd[2] / r1;
	npt.p[1] += V0 * sin(M_PI * r1 / g.LLn) * (crd[1] + .5 * g.LLy) / r1;
      }
    }
    if (r2 < g.LLn) {
      npt.n += (g.nn0 - g.nnb) * sqr(cos(M_PI / 2. * r2 / g.LLn));
      if (r2 > 0.0) {
	npt.p[2] += V0 * sin(M_PI * r2 / g.LLn) * crd[2] / r2;
	npt.p[1] += V0 * sin(M_PI * r2 / g.LLn) * (crd[1] - .5 * g.LLy) / r2;
      }
    }
    
    switch (kind) {
    case 0: // electrons
      // electron drift consistent with initial current
      if ((r1 <= g.LLn) && (r1 >= g.LLn - 2.*g.LLB)) {
	npt.p[0] = - g.BB * M_PI/(2.*g.LLB) * cos(M_PI * (g.LLn-r1)/(2.*g.LLB)) / npt.n;
      }
      if ((r2 <= g.LLn) && (r2 >= g.LLn - 2.*g.LLB)) {
	npt.p[0] = - g.BB * M_PI/(2.*g.LLB) * cos(M_PI * (g.LLn-r2)/(2.*g.LLB)) / npt.n;
      }
      
      npt.T[0] = g.TTe;
      npt.T[1] = g.TTe;
      npt.T[2] = g.TTe;
      break;
    case 1: // ions
      npt.T[0] = g.TTi;
      npt.T[1] = g.TTi;
      npt.T[2] = g.TTi;
      break;
    default:
      assert(0);
    }
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_partition
  
  std::vector<uint> setup_initial_partition()
  {
    SetupParticles<Mparticles> setup_particles;    
    return setup_particles.setup_partition(grid(), [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_particles
  
  void setup_initial_particles(Mparticles& mprts, std::vector<uint>& n_prts_by_patch)
  {
    SetupParticles<Mparticles> setup_particles;
    setup_particles.setup_particles(mprts, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(MfieldsState& mflds)
  {
    setupFields(grid(), mflds, [&](int m, double crd[3]) {
	double z1 = crd[2];
	double y1 = crd[1] + .5 * g.LLy;
	double r1 = sqrt(sqr(z1) + sqr(y1));
	double z2 = crd[2];
	double y2 = crd[1] - .5 * g.LLy;
	double r2 = sqrt(sqr(z2) + sqr(y2));

	double rv = 0.;
	switch (m) {
	case HZ:
	  if ( (r1 < g.LLn) && (r1 > g.LLn - 2*g.LLB) ) {
	    rv += - g.BB * sin(M_PI * (g.LLn - r1)/(2.*g.LLB)) * y1 / r1;
	  }
	  if ( (r2 < g.LLn) && (r2 > g.LLn - 2*g.LLB) ) {
	    rv += - g.BB * sin(M_PI * (g.LLn - r2)/(2.*g.LLB)) * y2 / r2;
	  }
	  return rv;
	  
	case HY:
	  if ( (r1 < g.LLn) && (r1 > g.LLn - 2*g.LLB) ) {
	    rv += g.BB * sin(M_PI * (g.LLn - r1)/(2.*g.LLB)) * z1 / r1;
	  }
	  if ( (r2 < g.LLn) && (r2 > g.LLn - 2*g.LLB) ) {
	    rv += g.BB * sin(M_PI * (g.LLn - r2)/(2.*g.LLB)) * z2 / r2;
	  }
	  return rv;
	  
	case EX:
	  if ( (r1 < g.LLn) && (r1 > g.LLn - 2*g.LLB) ) {
	    rv += g.MMach * sqrt(g.TTe/g.MMi) * g.BB *
	      sin(M_PI * (g.LLn - r1)/(2.*g.LLB)) * sin(M_PI * r1 / g.LLn);
	  }
	  if ( (r2 < g.LLn) && (r2 > g.LLn - 2*g.LLB) ) {
	    rv += g.MMach * sqrt(g.TTe/g.MMi) * g.BB *
	      sin(M_PI * (g.LLn - r2)/(2.*g.LLB)) * sin(M_PI * r2 / g.LLn);
	  }
	  return rv;

	  // FIXME, JXI isn't really needed anymore (?)
	case JXI:
	  if ( (r1 < g.LLn) && (r1 > g.LLn - 2*g.LLB) ) {
	    rv += g.BB * M_PI/(2.*g.LLB) * cos(M_PI * (g.LLn - r1)/(2.*g.LLB));
	  }
	  if ( (r2 < g.LLn) && (r2 > g.LLn - 2*g.LLB) ) {
	    rv += g.BB * M_PI/(2.*g.LLB) * cos(M_PI * (g.LLn - r2)/(2.*g.LLB));
	  }
	  return rv;
	  
	default:
	  return 0.;
	}
      });
  }
};

// ======================================================================
// run

static void run()
{
  g.BB = .07;
  g.nnb = .1;
  g.nn0 = 1.;
  g.MMach = 3.;
  g.LLn = 200.;
  g.LLB = 200./6.;
  g.TTe = .02;
  g.TTi = .02;
  g.MMi = 100.;
  
  g.LLy = 2. * g.LLn;
  g.LLz = 3. * g.LLn;

  psc_params.nmax = 1000; //32000;
  psc_params.stats_every = 100;
  
  
  PscBubble psc{psc_params};

  psc.initialize();
  psc.integrate();
}

// ======================================================================
// main

int
main(int argc, char **argv)
{
  psc_init(argc, argv);

  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
