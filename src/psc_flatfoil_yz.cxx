
#define DO_CUDA 1

#include <psc.h>

#ifdef USE_VPIC
#include "../libpsc/vpic/vpic_iface.h"
#endif

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
#include <marder.hxx>
#include <inject.hxx>
#include <heating.hxx>
#include <setup_particles.hxx>
#include <setup_fields.hxx>

#include "psc_particles_double.h"
#include "psc_fields_c.h"

#include "../libpsc/psc_sort/psc_sort_impl.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"
#include "../libpsc/psc_push_particles/push_config.hxx"
#include "../libpsc/psc_push_particles/push_dispatch.hxx"
#include "../libpsc/psc_push_particles/1vb/push_particles_1vbec_single.hxx"
#include "psc_push_fields_impl.hxx"
#include "bnd_particles_impl.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "../libpsc/psc_bnd_fields/psc_bnd_fields_impl.hxx"
#include "../libpsc/psc_inject/psc_inject_impl.hxx"
#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "../libpsc/psc_balance/psc_balance_impl.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"
#include "../libpsc/psc_push_fields/marder_impl.hxx"

#ifdef DO_CUDA
#include "../libpsc/cuda/push_particles_cuda_impl.hxx"
#include "../libpsc/cuda/push_fields_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_cuda_impl.hxx"
#include "../libpsc/cuda/bnd_particles_cuda_impl.hxx"
#include "../libpsc/cuda/inject_cuda_impl.hxx"
#include "../libpsc/cuda/heating_cuda_impl.hxx"
#include "../libpsc/cuda/checks_cuda_impl.hxx"
#include "../libpsc/cuda/marder_cuda_impl.hxx"
#include "../libpsc/cuda/sort_cuda_impl.hxx"
#include "../libpsc/cuda/collision_cuda_impl.hxx"
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
#endif

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// ======================================================================
// InjectFoil

struct InjectFoilParams
{
  double yl, yh;
  double zl, zh;
  double n;
  double Te, Ti;
};

struct InjectFoil : InjectFoilParams
{
  InjectFoil() = default;
  
  InjectFoil(const InjectFoilParams& params)
    : InjectFoilParams(params)
  {}

  bool is_inside(double crd[3])
  {
    return (crd[1] >= yl && crd[1] <= yh &&
	    crd[2] >= zl && crd[2] <= zh);
  }

  void init_npt(int pop, double crd[3], struct psc_particle_npt *npt)
  {
    if (!is_inside(crd)) {
      npt->n = 0;
      return;
    }
    
    switch (pop) {
    case MY_ION:
      npt->n    = n;
      npt->T[0] = Ti;
      npt->T[1] = Ti;
      npt->T[2] = Ti;
      break;
    case MY_ELECTRON:
      npt->n    = n;
      npt->T[0] = Te;
      npt->T[1] = Te;
      npt->T[2] = Te;
      break;
    default:
      assert(0);
    }
  }
};

// ======================================================================
// HeatingSpotFoil

struct HeatingSpotFoilParams
{
  double zl; // in internal units (d_e)
  double zh;
  double xc;
  double yc;
  double rH;
  double T;
  double Mi;
};

struct HeatingSpotFoil : HeatingSpotFoilParams
{
  HeatingSpotFoil() = default;
  
  HeatingSpotFoil(const HeatingSpotFoilParams& params)
    : HeatingSpotFoilParams(params)
  {
    double width = zh - zl;
    fac = (8.f * pow(T, 1.5)) / (sqrt(Mi) * width);
    // FIXME, I don't understand the sqrt(Mi) in here
  }
  
  double operator()(const double *crd)
  {
    double x = crd[0], y = crd[1], z = crd[2];

    if (z <= zl || z >= zh) {
      return 0;
    }
    
    return fac * exp(-(sqr(x-xc) + sqr(y-yc)) / sqr(rH));
  }

private:
  double fac;
};

// ======================================================================
// PscFlatfoilParams

struct PscFlatfoilParams
{
  double BB;
  double Zi;

  double background_n;
  double background_Te;
  double background_Ti;

  int sort_interval;

  int collision_interval;
  double collision_nu;

  int marder_interval;
  double marder_diffusion;
  int marder_loop;
  bool marder_dump;

  int balance_interval;
  double balance_factor_fields;
  bool balance_print_loads;
  bool balance_write_loads;

  bool inject_enable;
  int inject_kind_n;
  int inject_interval;
  int inject_tau;
  InjectFoil inject_target;

  int heating_begin;
  int heating_end;
  int heating_interval;

  ChecksParams checks_params;
};

// ======================================================================
// PscFlatfoil
//
// eventually, a Psc replacement / derived class, but for now just
// pretending to be something like that

struct PscFlatfoil : PscFlatfoilParams
{
#if DO_CUDA
  using Mparticles_t = MparticlesCuda;
  using Mfields_t = MfieldsCuda;
  using Config1vbec3d = PushParticlesConfig<BS144, IpEc, DepositVb3d, CurrmemShared>;
  using PushParticles_t = PushParticlesCuda<Config1vbec3d>;
  using Sort_t = SortCuda;
  using Collision_t = CollisionCuda;
  using PushFields_t = PushFieldsCuda;
  using BndParticles_t = BndParticlesCuda;
  using Bnd_t = BndCuda;
  using BndFields_t = BndFieldsNone<Mfields_t>;
  using Inject_t = InjectCuda<InjectFoil>;
  using Heating_t = HeatingCuda;
  using Balance_t = Balance_<MparticlesSingle, MfieldsSingle>;
  using Checks_t = ChecksCuda;
  using Marder_t = MarderCuda;
#else
  using Mparticles_t = MparticlesDouble;
  using Mfields_t = MfieldsC;
#if 0 // generic_c: 2nd order
  using PushParticles_t = PushParticles__<Config2nd<dim_yz>>;
  using checks_order = checks_order_2nd;
#else // 1vbec: 1st order Villasenor-Buneman energy-conserving (kinda...)
  using PushParticles_t = PushParticles1vb<Config1vbec<Mparticles_t, Mfields_t, dim_yz>>;
  using checks_order = checks_order_1st;
#endif
  using Sort_t = SortCountsort2<Mparticles_t>;
  using Collision_t = Collision_<Mparticles_t, Mfields_t>;
  using PushFields_t = PushFields<Mfields_t>;
  using BndParticles_t = BndParticles_<Mparticles_t>;
  using Bnd_t = Bnd_<Mfields_t>;
  using BndFields_t = BndFieldsNone<Mfields_t>; // FIXME, why MfieldsC hardcoded???
  using Inject_t = Inject_<Mparticles_t, MfieldsC, InjectFoil>; // FIXME, shouldn't always use MfieldsC
  using Heating_t = Heating__<Mparticles_t>;
  using Balance_t = Balance_<Mparticles_t, Mfields_t>;
  using Checks_t = Checks_<Mparticles_t, Mfields_t, checks_order>;
  using Marder_t = Marder_<Mparticles_t, Mfields_t>;
#endif
  
  PscFlatfoil(const PscFlatfoilParams& params, Heating_t heating, psc *psc)
    : PscFlatfoilParams(params),
      psc_{psc},
      mprts_{dynamic_cast<Mparticles_t&>(*PscMparticlesBase{psc->particles}.sub())},
      mflds_{dynamic_cast<Mfields_t&>(*PscMfieldsBase{psc->flds}.sub())},
      collision_{psc_comm(psc), collision_interval, collision_nu},
      bndp_{psc_->mrc_domain_, psc_->grid()},
      bnd_{psc_->grid(), psc_->mrc_domain_, psc_->ibn},
      balance_{balance_interval, balance_factor_fields, balance_print_loads, balance_write_loads},
      heating_{heating},
      inject_{psc_comm(psc), inject_enable, inject_interval, inject_tau, inject_kind_n, inject_target},
      checks_{psc_comm(psc), checks_params},
      marder_(psc_comm(psc), marder_diffusion, marder_loop, marder_dump)
  {
    MPI_Comm comm = psc_comm(psc_);

    // --- partition particles and initial balancing
    mpi_printf(comm, "**** Partitioning...\n");
    auto n_prts_by_patch_old = setup_initial_partition();
    auto n_prts_by_patch_new = balance_.initial(psc_, n_prts_by_patch_old);
    // balance::initial does not rebalance particles, because the old way of doing this
    // does't even have the particle data structure created yet -- FIXME?
    mprts_.reset(psc_->grid());
    
    mpi_printf(comm, "**** Setting up particles...\n");
    setup_initial_particles(mprts_, n_prts_by_patch_new);
    
    mpi_printf(comm, "**** Setting up fields...\n");
    setup_initial_fields(mflds_);

    checks_.gauss(mprts_, mflds_);
    psc_setup_member_objs(psc_);

    setup_stats();
  }

  void init_npt(int kind, double crd[3], psc_particle_npt& npt)
  {
    switch (kind) {
    case MY_ION:
      npt.n    = background_n;
      npt.T[0] = background_Ti;
      npt.T[1] = background_Ti;
      npt.T[2] = background_Ti;
      break;
    case MY_ELECTRON:
      npt.n    = background_n;
      npt.T[0] = background_Te;
      npt.T[1] = background_Te;
      npt.T[2] = background_Te;
      break;
    default:
      assert(0);
    }
      
    if (inject_target.is_inside(crd)) {
      // replace values above by target values
      inject_target.init_npt(kind, crd, &npt);
    }
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_partition
  
  std::vector<uint> setup_initial_partition()
  {
    return SetupParticles<Mparticles_t>::setup_partition(psc_, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }
  
  // ----------------------------------------------------------------------
  // setup_initial_particles
  
  void setup_initial_particles(Mparticles_t& mprts, std::vector<uint>& n_prts_by_patch)
  {
    SetupParticles<Mparticles_t>::setup_particles(mprts, psc_, n_prts_by_patch, [&](int kind, double crd[3], psc_particle_npt& npt) {
	this->init_npt(kind, crd, npt);
      });
  }

  // ----------------------------------------------------------------------
  // setup_initial_fields
  
  void setup_initial_fields(Mfields_t& mflds)
  {
    SetupFields<Mfields_t>::set(mflds, [&](int m, double crd[3]) {
	switch (m) {
	case HY: return BB;
	default: return 0.;
	}
      });
  }

  // ----------------------------------------------------------------------
  // setup_stats
  
  void setup_stats()
  {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");

    // generic stats categories
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_comm = psc_stats_register("time communication");
    st_time_output = psc_stats_register("time output");
  }
  
  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    //psc_method_initialize(psc_->method, psc_);
    psc_output(psc_);
    psc_stats_log(psc_);
    psc_print_profiling(psc_);

    mpi_printf(psc_comm(psc_), "Initialization complete.\n");

    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(psc_comm(psc_), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (psc_->timestep < psc_->prm.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  psc_->prm.write_checkpoint_every_step > 0 &&
	  psc_->timestep % psc_->prm.write_checkpoint_every_step == 0) {
	psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(psc_comm(psc_), "**** Step %d / %d, Time %g\n", psc_->timestep + 1,
		 psc_->prm.nmax, psc_->timestep * psc_->dt);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
    
      psc_->timestep++; // FIXME, too hacky
      psc_output(psc_);

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_.get_n_prts();

      if (psc_->timestep % psc_->prm.stats_every == 0) {
	psc_stats_log(psc_);
	psc_print_profiling(psc_);
      }

      if (psc_->prm.wallclock_limit > 0.) {
	double wallclock_elapsed = MPI_Wtime() - psc_->time_start;
	double wallclock_elapsed_max;
	MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		      MPI_COMM_WORLD);
      
	if (wallclock_elapsed_max > psc_->prm.wallclock_limit) {
	  mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	  break;
	}
      }
    }

    if (psc_->prm.write_checkpoint) {
      psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(psc_comm(psc_), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

  // ----------------------------------------------------------------------
  // step
  //
  // things are missing from the generic step():
  // - pushp prep

  void step()
  {
    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject, pr_heating;
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
    }

    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = psc_comm(psc_);
    int timestep = psc_->timestep;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      balance_(psc_, mprts_);
    }

    if (sort_interval > 0 && timestep % sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      sort_(mprts_);
      prof_stop(pr_sort);
    }
    
    if (collision_interval > 0 && timestep % collision_interval == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      collision_(mprts_);
      prof_stop(pr_collision);
    }
    
    if (checks_params.continuity_every_step > 0 && timestep % checks_params.continuity_every_step == 0) {
      prof_start(pr_checks);
      checks_.continuity_before_particle_push(mprts_);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_.push_mprts(mprts_, mflds_);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}
    
    prof_start(pr_bndp);
    bndp_(mprts_);
    prof_stop(pr_bndp);
    
    // === field propagation B^{n+1/2} -> B^{n+1}
    prof_start(pr_push_flds);
    pushf_.push_H(mflds_, .5, dim_yz{});
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}
    
    prof_start(pr_inject);
    inject_(mprts_);
    prof_stop(pr_inject);
      
    // only heating between heating_tb and heating_te
    if (timestep >= heating_begin && timestep < heating_end &&
	timestep % heating_interval == 0) {
      prof_start(pr_heating);
      heating_(mprts_);
      prof_stop(pr_heating);
    }

    // === field propagation E^{n+1/2} -> E^{n+3/2}
    prof_start(pr_bndf);
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    
    bndf_.add_ghosts_J(mflds_);
    bnd_.add_ghosts(mflds_, JXI, JXI + 3);
    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);
    prof_stop(pr_bndf);
    
    prof_restart(pr_push_flds);
    pushf_.push_E(mflds_, 1., dim_yz{});
    prof_stop(pr_push_flds);
    
    prof_restart(pr_bndf);
    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
      
    // === field propagation B^{n+1} -> B^{n+3/2}
    prof_restart(pr_push_flds);
    pushf_.push_H(mflds_, .5, dim_yz{});
    prof_stop(pr_push_flds);
    
    prof_start(pr_bndf);
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
      
    if (checks_params.continuity_every_step > 0 && timestep % checks_params.continuity_every_step == 0) {
      prof_restart(pr_checks);
      checks_.continuity_after_particle_push(mprts_, mflds_);
      prof_stop(pr_checks);
    }
    
    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      marder_(mflds_, mprts_);
      prof_stop(pr_marder);
    }
    
    if (checks_params.gauss_every_step > 0 && timestep % checks_params.gauss_every_step == 0) {
      prof_restart(pr_checks);
      checks_.gauss(mprts_, mflds_);
      prof_stop(pr_checks);
    }
    
    //psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

private:
  psc* psc_;
  Mparticles_t& mprts_;
  Mfields_t& mflds_;

  Sort_t sort_;
  Collision_t collision_;
  PushParticles_t pushp_;
  PushFields_t pushf_;
  BndParticles_t bndp_;
  Bnd_t bnd_;
  BndFields_t bndf_;
  Balance_t balance_;

  Heating_t heating_;
  Inject_t inject_;

  Checks_t checks_;
  Marder_t marder_;
  
  int st_nr_particles;
  int st_time_step;
};

// ======================================================================
// PscFlatfoilBuilder

struct PscFlatfoilBuilder
{
  using Heating_t = PscFlatfoil::Heating_t;

  PscFlatfoilBuilder()
    : psc_(psc_create(MPI_COMM_WORLD))
  {}
  
  PscFlatfoil* makePscFlatfoil();

  PscFlatfoilParams params;
  
  // state
  double d_i;
  double LLs;
  double LLn;

  psc* psc_;
};

// ----------------------------------------------------------------------
// PscFlatfoilBuilder::makePscFlatfoil

PscFlatfoil* PscFlatfoilBuilder::makePscFlatfoil()
{
  MPI_Comm comm = psc_comm(psc_);
  
  mpi_printf(comm, "*** Setting up...\n");

  psc_default_dimensionless(psc_);

  psc_->prm.nmax = 201;
  psc_->prm.nicell = 10;
  psc_->prm.fractional_n_particles_per_cell = true;
  psc_->prm.cfl = 0.75;

  // --- setup domain
  double LLy = 400.;
  double LLz = 400.*4.;

  auto grid_domain = Grid_t::Domain{{1, 32, 128}, // global number of grid points
				    {1., LLy, LLz}, {0., -.5*LLy, -.5*LLz}, // domain size, origin
				    {1, 1, 4}, // division into patches
                                    {1, 4, 4}}; // blocksize FIXME needs to match cuda pusher

  auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			{ BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			{ BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

  psc_set_from_options(psc_);

  params.BB = 0.;
  params.Zi = 1.;

  // --- for background plasma
  params.background_n  = .002;
  params.background_Te = .001;
  params.background_Ti = .001;
  
  // -- setup particles
  // last population ("e") is neutralizing
 // FIXME, hardcoded mass ratio 100
  Grid_t::Kinds kinds = {{params.Zi, 100.*params.Zi, "i"}, { -1., 1., "e"}};
  psc_->prm.neutralizing_population = MY_ELECTRON;
  
  d_i = sqrt(kinds[MY_ION].m / kinds[MY_ION].q);

  mpi_printf(comm, "d_e = %g, d_i = %g\n", 1., d_i);
  mpi_printf(comm, "lambda_De (background) = %g\n", sqrt(params.background_Te));

  // sort
  params.sort_interval = 10;

  // collisions
  params.collision_interval = 10;
  params.collision_nu = .1;

  // --- setup heating
  double heating_zl = -1.;
  double heating_zh =  1.;
  double heating_xc = 0.;
  double heating_yc = 0.;
  double heating_rH = 3.;
  auto heating_foil_params = HeatingSpotFoilParams{};
  heating_foil_params.zl = heating_zl * d_i;
  heating_foil_params.zh = heating_zh * d_i;
  heating_foil_params.xc = heating_xc * d_i;
  heating_foil_params.yc = heating_yc * d_i;
  heating_foil_params.rH = heating_rH * d_i;
  heating_foil_params.T  = .04;
  heating_foil_params.Mi = heating_rH * kinds[MY_ION].m;
  auto heating_spot = HeatingSpotFoil{heating_foil_params};
  params.heating_interval = 20;
  params.heating_begin = 0;
  params.heating_end = 10000000;
  int heating_kind = MY_ELECTRON;
  auto heating = Heating_t{params.heating_interval, heating_kind, heating_spot};

  // -- setup injection
  double target_yl     = -100000.;
  double target_yh     =  100000.;
  double target_zwidth =  1.;
  auto inject_foil_params = InjectFoilParams{};
  inject_foil_params.yl =   target_yl * d_i;
  inject_foil_params.yh =   target_yh * d_i;
  inject_foil_params.zl = - target_zwidth * d_i;
  inject_foil_params.zh =   target_zwidth * d_i;
  inject_foil_params.n  = 1.;
  inject_foil_params.Te = .001;
  inject_foil_params.Ti = .001;
  params.inject_target = InjectFoil{inject_foil_params};
  params.inject_enable = true;
  params.inject_kind_n = MY_ELECTRON;
  params.inject_interval = 20;
  params.inject_tau = 40;

  // --- checks
  params.checks_params.continuity_every_step = -1;
  params.checks_params.continuity_threshold = 1e-14;
  params.checks_params.continuity_verbose = true;
  params.checks_params.continuity_dump_always = false;

  params.checks_params.gauss_every_step = -1;
  params.checks_params.gauss_threshold = 1e-14;
  params.checks_params.gauss_verbose = true;
  params.checks_params.gauss_dump_always = false;

  // --- marder
  params.marder_interval = 5;
  params.marder_diffusion = 0.9;
  params.marder_loop = 3;
  params.marder_dump = false;

  // --- balancing
  params.balance_interval = 100;
  params.balance_factor_fields = 0.1;
  params.balance_print_loads = true;
  params.balance_write_loads = false;

  // --- generic setup
  psc_setup_coeff(psc_);
  psc_setup_domain(psc_, grid_domain, grid_bc, kinds);

  // make sure that np isn't overridden on the command line
  Int3 np;
  mrc_domain_get_param_int3(psc_->mrc_domain_, "np", np);
  assert(np == grid_domain.np);

  // --- create and initialize base particle data structure x^{n+1/2}, p^{n+1/2}
  mpi_printf(comm, "**** Creating particle data structure...\n");
  psc_->particles = PscMparticlesCreate(comm, psc_->grid(),
					Mparticles_traits<PscFlatfoil::Mparticles_t>::name).mprts();

  // --- create and set up base mflds
  mpi_printf(comm, "**** Creating fields...\n");
  psc_->flds = PscMfieldsCreate(comm, psc_->grid(), psc_->n_state_fields, psc_->ibn,
				Mfields_traits<PscFlatfoil::Mfields_t>::name).mflds();

  return new PscFlatfoil(params, heating, psc_);
}

// ----------------------------------------------------------------------
// psc_ops "flatfoil"

struct psc_ops_flatfoil : psc_ops {
  psc_ops_flatfoil() {
    name             = "flatfoil";
  }
} psc_flatfoil_ops;

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

  mrc_class_register_subclass(&mrc_class_psc, &psc_flatfoil_ops);

  auto sim = PscFlatfoilBuilder{};
  auto flatfoil = sim.makePscFlatfoil();
  
  psc_view(sim.psc_);
  psc_mparticles_view(sim.psc_->particles);
  psc_mfields_view(sim.psc_->flds);
  
  flatfoil->integrate();
  
  delete flatfoil;
  
  psc_destroy(sim.psc_);
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}

