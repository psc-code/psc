
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
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

#include "psc_config.hxx"

#include <mrc_params.h>

#include <math.h>

struct psc_bubble {
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

#define to_psc_bubble(psc) mrc_to_subobj(psc, struct psc_bubble)

#define VAR(x) (void *)offsetof(struct psc_bubble, x)
static struct param psc_bubble_descr[] = {
  { "BB"            , VAR(BB)              , PARAM_DOUBLE(.07)    },
  { "nnb"           , VAR(nnb)             , PARAM_DOUBLE(.1)     },
  { "nn0"           , VAR(nn0)             , PARAM_DOUBLE(1.)     },
  { "MMach"         , VAR(MMach)           , PARAM_DOUBLE(3.)     },
  { "LLn"           , VAR(LLn)             , PARAM_DOUBLE(200.)   },
  { "LLB"           , VAR(LLB)             , PARAM_DOUBLE(200./6.)},
  { "LLz"           , VAR(LLz)             , PARAM_DOUBLE(0.)     },
  { "LLy"           , VAR(LLy)             , PARAM_DOUBLE(0.)     },
  { "TTe"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "TTi"           , VAR(TTe)             , PARAM_DOUBLE(.02)    },
  { "MMi"           , VAR(MMi)             , PARAM_DOUBLE(100.)   },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_bubble_create

static void
psc_bubble_create(struct psc *psc)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);
  psc_default_dimensionless(psc);

  psc->prm.nmax = 1000; //32000;
  psc->prm.nicell = 100;

  bubble->LLy = 2. * bubble->LLn;
  bubble->LLz = 3. * bubble->LLn;

  psc->domain_ = Grid_t::Domain{{1, 128, 512},
				{bubble->LLn, bubble->LLy, bubble->LLz},
				{0., -.5 * bubble->LLy, -.5 * bubble->LLz},
				{1, 1, 4}};
  
  psc->bc_ = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
		    { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
		    { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
		    { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");
}

// ----------------------------------------------------------------------
// psc_bubble_setup

static void
psc_bubble_setup(struct psc *psc)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);
  
  psc_setup_super(psc);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "lambda_D = %g\n", sqrt(bubble->TTe));
}

// ----------------------------------------------------------------------
// psc_bubble_read

static void
psc_bubble_read(struct psc *psc, struct mrc_io *io)
{
  psc_read_super(psc, io);
}

// ----------------------------------------------------------------------
// psc_bubble_init_field

static double
psc_bubble_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  double BB = bubble->BB;
  double LLn = bubble->LLn;
  double LLy = bubble->LLy;
  double LLB = bubble->LLB;
  double MMi = bubble->MMi;
  double MMach = bubble->MMach;
  double TTe = bubble->TTe;

  double z1 = x[2];
  double y1 = x[1] + .5 * LLy;
  double r1 = sqrt(sqr(z1) + sqr(y1));
  double z2 = x[2];
  double y2 = x[1] - .5 * LLy;
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
}

// ----------------------------------------------------------------------
// psc_bubble_init_npt

static void
psc_bubble_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  struct psc_bubble *bubble = to_psc_bubble(psc);

  double BB = bubble->BB;
  double LLy = bubble->LLy;
  double LLn = bubble->LLn;
  double LLB = bubble->LLB;
  double V0 = bubble->MMach * sqrt(bubble->TTe / bubble->MMi);

  double nnb = bubble->nnb;
  double nn0 = bubble->nn0;

  double TTe = bubble->TTe, TTi = bubble->TTi;

  double r1 = sqrt(sqr(x[2]) + sqr(x[1] + .5 * LLy));
  double r2 = sqrt(sqr(x[2]) + sqr(x[1] - .5 * LLy));

  npt->n = nnb;
  if (r1 < LLn) {
    npt->n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r1 / LLn));
    if (r1 > 0.0) {
      npt->p[2] += V0 * sin(M_PI * r1 / LLn) * x[2] / r1;
      npt->p[1] += V0 * sin(M_PI * r1 / LLn) * (x[1] + .5 * LLy) / r1;
    }
  }
  if (r2 < LLn) {
    npt->n += (nn0 - nnb) * sqr(cos(M_PI / 2. * r2 / LLn));
    if (r2 > 0.0) {
      npt->p[2] += V0 * sin(M_PI * r2 / LLn) * x[2] / r2;
      npt->p[1] += V0 * sin(M_PI * r2 / LLn) * (x[1] - .5 * LLy) / r2;
    }
  }

  switch (kind) {
  case 0: // electrons
    // electron drift consistent with initial current
    if ((r1 <= LLn) && (r1 >= LLn - 2.*LLB)) {
      npt->p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r1)/(2.*LLB)) / npt->n;
    }
    if ((r2 <= LLn) && (r2 >= LLn - 2.*LLB)) {
      npt->p[0] = - BB * M_PI/(2.*LLB) * cos(M_PI * (LLn-r2)/(2.*LLB)) / npt->n;
    }

    npt->T[0] = TTe;
    npt->T[1] = TTe;
    npt->T[2] = TTe;
    break;
  case 1: // ions
    npt->T[0] = TTi;
    npt->T[1] = TTi;
    npt->T[2] = TTi;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_bubble_ops

struct psc_ops_bubble : psc_ops {
  psc_ops_bubble() {
    name             = "bubble";
    size             = sizeof(struct psc_bubble);
    param_descr      = psc_bubble_descr;
    create           = psc_bubble_create;
    setup            = psc_bubble_setup;
    read             = psc_bubble_read;
    init_field       = psc_bubble_init_field;
    init_npt         = psc_bubble_init_npt;
  }
} psc_bubble_ops;

// ======================================================================
// PscBubbleParams

struct PscBubbleParams
{
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

  ChecksParams checks_params;
};

using PscConfig = PscConfig1vbecSingle<dim_yz>;

// ======================================================================
// PscBubble

struct PscBubble : PscBubbleParams
{
  using Mparticles_t = PscConfig::Mparticles_t;
  using Mfields_t = PscConfig::Mfields_t;
  using Sort_t = PscConfig::Sort_t;
  using Collision_t = PscConfig::Collision_t;
  using PushParticles_t = PscConfig::PushParticles_t;
  using PushFields_t = PscConfig::PushFields_t;

  PscBubble(const PscBubbleParams& params, psc *psc)
    : PscBubbleParams(params),
      psc_{psc},
      mprts_{dynamic_cast<Mparticles_t&>(*PscMparticlesBase{psc->particles}.sub())},
      mflds_{dynamic_cast<Mfields_t&>(*PscMfieldsBase{psc->flds}.sub())},
      collision_{psc_comm(psc), collision_interval, collision_nu}
  {
    setup_stats();
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

      mpi_printf(psc_comm(psc_), "**** Step %d / %d, Code Time %g, Wall Time %g\n", psc_->timestep + 1,
		 psc_->prm.nmax, psc_->timestep * psc_->dt, MPI_Wtime() - psc_->time_start);

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
    PscMparticlesBase mprts(psc_->particles);
    PscMfieldsBase mflds(psc_->flds);
    PscPushFieldsBase pushf(psc_->push_fields);
    PscBndParticlesBase bndp(psc_->bnd_particles);
    PscBndBase bnd(psc_->bnd);
    PscBndFieldsBase bndf(pushf.pushf()->bnd_fields);
    auto balance = PscBalanceBase{psc_->balance};

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
    MPI_Comm comm = psc_comm(psc_);
    int timestep = psc_->timestep;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      balance(psc_, mprts);
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
      mpi_printf(comm, "***** Checking continuity...\n");
      prof_start(pr_checks);
      PscChecksBase{psc_->checks}.continuity_before_particle_push(psc_);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_.push_mprts(mprts_, mflds_);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

#if 0
    prof_start(pr_sync1);
    MPI_Barrier(comm);
    prof_stop(pr_sync1);
#endif
    
    // === field propagation B^{n+1/2} -> B^{n+1}
    prof_start(pr_push_flds);
    pushf_.push_H(mflds, .5);
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_sync3);
    MPI_Barrier(comm);
    prof_stop(pr_sync3);

    prof_start(pr_bndp);
    bndp(*mprts.sub());
    prof_stop(pr_bndp);

    // === field propagation E^{n+1/2} -> E^{n+3/2}
#if 1
    prof_start(pr_bndf);
    bndf.fill_ghosts_H(mflds);
    bnd.fill_ghosts(mflds, HX, HX + 3);
#endif

    bndf.add_ghosts_J(mflds);
    bnd.add_ghosts(mflds, JXI, JXI + 3);
    bnd.fill_ghosts(mflds, JXI, JXI + 3);
    prof_stop(pr_bndf);
    
#if 1
    prof_start(pr_sync4a);
    MPI_Barrier(comm);
    prof_stop(pr_sync4a);
#endif
    
    prof_restart(pr_push_flds);
    pushf_.push_E(mflds, 1.);
    prof_stop(pr_push_flds);
    
#if 0
    prof_start(pr_sync4b);
    MPI_Barrier(comm);
    prof_stop(pr_sync4b);
#endif

#if 1
    prof_restart(pr_bndf);
    bndf.fill_ghosts_E(mflds);
    bnd.fill_ghosts(mflds, EX, EX + 3);
    prof_stop(pr_bndf);
#endif
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
      
    // === field propagation B^{n+1} -> B^{n+3/2}
    prof_restart(pr_push_flds);
    pushf_.push_H(mflds, .5);
    prof_stop(pr_push_flds);

#if 1
    prof_start(pr_bndf);
    bndf.fill_ghosts_H(mflds);
    bnd.fill_ghosts(mflds, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
#endif

#if 0
    prof_start(pr_sync5);
    MPI_Barrier(comm);
    prof_stop(pr_sync5);
#endif
    
    if (checks_params.continuity_every_step > 0 && timestep % checks_params.continuity_every_step == 0) {
      prof_restart(pr_checks);
      PscChecksBase{psc_->checks}.continuity_after_particle_push(psc_);
      prof_stop(pr_checks);
    }
    
    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      PscMarderBase{psc_->marder}(mflds, mprts);
      prof_stop(pr_marder);
    }
    
    if (checks_params.gauss_every_step > 0 && timestep % checks_params.gauss_every_step == 0) {
      prof_restart(pr_checks);
      PscChecksBase{psc_->checks}.gauss(psc_);
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
  
  int st_nr_particles;
  int st_time_step;
};

// ======================================================================
// PscBubbleBuilder

struct PscBubbleBuilder
{
  PscBubbleBuilder()
    : psc_(psc_create(MPI_COMM_WORLD))
  {}

  PscBubble* makePscBubble();

  PscBubbleParams params;
  psc* psc_;
};

// ----------------------------------------------------------------------
// PscBubbleBuilder::makePscBubble

PscBubble* PscBubbleBuilder::makePscBubble()
{
  MPI_Comm comm = psc_comm(psc_);
  
  mpi_printf(comm, "*** Setting up...\n");

  // sort
  params.sort_interval = 10;

  // collisions
  params.collision_interval = 10;
  params.collision_nu = .1;

  // --- checks
  params.checks_params.continuity_every_step = -1;
  params.checks_params.continuity_threshold = 1e-6;
  params.checks_params.continuity_verbose = true;
  params.checks_params.continuity_dump_always = false;

  params.checks_params.gauss_every_step = -1;
  params.checks_params.gauss_threshold = 1e-6;
  params.checks_params.gauss_verbose = true;
  params.checks_params.gauss_dump_always = false;

  // --- marder
  params.marder_interval = 0*5;
  params.marder_diffusion = 0.9;
  params.marder_loop = 3;
  params.marder_dump = false;

  // --- balancing
  params.balance_interval = 0;
  params.balance_factor_fields = 0.1;
  params.balance_print_loads = true;
  params.balance_write_loads = false;
  
  psc_set_from_options(psc_);
  psc_setup(psc_);

  return new PscBubble{params, psc_};
}

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

  mrc_class_register_subclass(&mrc_class_psc, &psc_bubble_ops);

  auto sim = PscBubbleBuilder{};
  auto bubble = sim.makePscBubble();

  psc_view(sim.psc_);
  psc_mparticles_view(sim.psc_->particles);
  psc_mfields_view(sim.psc_->flds);
  
  bubble->integrate();

  delete bubble;
  
  psc_destroy(sim.psc_);
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
