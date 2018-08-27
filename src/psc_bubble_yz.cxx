
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
    double balance_factor_fields = .1;
    bool balance_print_loads = true;
    bool balance_write_loads = false;
    balance_.reset(new Balance_t{balance_interval, balance_factor_fields,
	  balance_print_loads, balance_write_loads});

    // -- Sort
    // FIXME, needs a way to make it gets set?
    sort_interval = 10;

    // -- Collision
    int collision_interval = 10;
    double collision_nu = .1;
    collision_.reset(new Collision_t{grid(), collision_interval, collision_nu});

    // -- Checks
    ChecksParams checks_params;

    checks_params.continuity_every_step = -1;
    checks_params.continuity_threshold = 1e-6;
    checks_params.continuity_verbose = true;
    checks_params.continuity_dump_always = false;
    
    checks_params.gauss_every_step = -1;
    checks_params.gauss_threshold = 1e-6;
    checks_params.gauss_verbose = true;
    checks_params.gauss_dump_always = false;

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
    auto n_prts_by_patch_new = balance_->initial(psc_, n_prts_by_patch_old);
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
    setup_particles.setup_particles(mprts, psc_, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
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

  // ----------------------------------------------------------------------
  // step
  //
  // things are missing from the generic step():
  // - pushp prep

  void step()
  {
    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject, pr_heating,
      pr_sync1, pr_sync2, pr_sync3, pr_sync4, pr_sync5, pr_sync4a, pr_sync4b;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
      pr_inject = prof_register("step_inject", 1., 0, 0);
      pr_heating = prof_register("step_heating", 1., 0, 0);
      pr_sync1 = prof_register("step_sync1", 1., 0, 0);
      pr_sync2 = prof_register("step_sync2", 1., 0, 0);
      pr_sync3 = prof_register("step_sync3", 1., 0, 0);
      pr_sync4 = prof_register("step_sync4", 1., 0, 0);
      pr_sync5 = prof_register("step_sync5", 1., 0, 0);
      pr_sync4a = prof_register("step_sync4a", 1., 0, 0);
      pr_sync4b = prof_register("step_sync4b", 1., 0, 0);
    }

    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = grid().comm();
    int timestep = grid().timestep();

    auto& mprts = *mprts_;
    auto& mflds = *mflds_;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      (*balance_)(psc_, mprts);
    }

    if (sort_interval > 0 && timestep % sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      (*sort_)(mprts);
      prof_stop(pr_sort);
    }
    
    if (collision_->interval() > 0 && timestep % collision_->interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      (*collision_)(mprts);
      prof_stop(pr_collision);
    }
    
    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      mpi_printf(comm, "***** Checking continuity...\n");
      prof_start(pr_checks);
      checks_->continuity_before_particle_push(mprts);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_->push_mprts(mprts, mflds);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

#if 0
    prof_start(pr_sync1);
    MPI_Barrier(comm);
    prof_stop(pr_sync1);
#endif
    
    // === field propagation B^{n+1/2} -> B^{n+1}
    prof_start(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_sync3);
    MPI_Barrier(comm);
    prof_stop(pr_sync3);

    prof_start(pr_bndp);
    (*bndp_)(mprts);
    prof_stop(pr_bndp);

    // === field propagation E^{n+1/2} -> E^{n+3/2}
#if 1
    prof_start(pr_bndf);
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
#endif

    bndf_->add_ghosts_J(mflds);
    bnd_->add_ghosts(mflds, JXI, JXI + 3);
    bnd_->fill_ghosts(mflds, JXI, JXI + 3);
    prof_stop(pr_bndf);
    
#if 1
    prof_start(pr_sync4a);
    MPI_Barrier(comm);
    prof_stop(pr_sync4a);
#endif
    
    prof_restart(pr_push_flds);
    pushf_->push_E(mflds, 1., DIM{});
    prof_stop(pr_push_flds);
    
#if 0
    prof_start(pr_sync4b);
    MPI_Barrier(comm);
    prof_stop(pr_sync4b);
#endif

#if 1
    prof_restart(pr_bndf);
    bndf_->fill_ghosts_E(mflds);
    bnd_->fill_ghosts(mflds, EX, EX + 3);
    prof_stop(pr_bndf);
#endif
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
      
    // === field propagation B^{n+1} -> B^{n+3/2}
    prof_restart(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);

#if 1
    prof_start(pr_bndf);
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
#endif

#if 0
    prof_start(pr_sync5);
    MPI_Barrier(comm);
    prof_stop(pr_sync5);
#endif
    
    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      prof_restart(pr_checks);
      checks_->continuity_after_particle_push(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      (*marder_)(mflds, mprts);
      prof_stop(pr_marder);
    }
    
    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      prof_restart(pr_checks);
      checks_->gauss(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    //psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
#ifdef USE_VPIC
  vpic_base_init(&argc, &argv);
#else
  MPI_Init(&argc, &argv);
#endif
  libmrc_params_init(argc, argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  auto psc = new PscBubble;

  psc->initialize();
  psc->integrate();

  delete psc;
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
