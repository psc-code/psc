
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_heating.h>
#include <psc_inject.h>
#include <psc_target_private.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define USE_OWN_PSC_STEP

#include <psc_balance.h>
#include <psc_sort.h>
#include <psc_collision.h>
#include <psc_checks.h>
#include <psc_bnd_particles.h>
#include <psc_marder.h>
#include <psc_method.h>

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
#include <inject.hxx>
#include <heating.hxx>

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

// ======================================================================
// psc subclass "flatfoil"

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

struct psc_flatfoil {
  double BB;
  double Zi;
  double LLf;
  double LLz;
  double LLy;

  double background_n;
  double background_Te;
  double background_Ti;

  bool no_initial_target; // for testing, the target can be turned off in the initial condition

  double target_yl;
  double target_yh;
  double target_zwidth;
  struct psc_target *target;

  struct psc_inject *inject;

  double heating_zl; // this is ugly as these are used to set the corresponding
  double heating_zh; // quantities in psc_heating, but having them here we can rescale
  double heating_xc; // them from d_i to internal (d_e) units
  double heating_yc;
  double heating_rH;
  
  // state
  double d_i;
  double LLs;
  double LLn;
};

#define psc_flatfoil(psc) mrc_to_subobj(psc, struct psc_flatfoil)

#define VAR(x) (void *)offsetof(struct psc_flatfoil, x)
static struct param psc_flatfoil_descr[] = {
  { "BB"                , VAR(BB)                , PARAM_DOUBLE(.0)         },
  { "Zi"                , VAR(Zi)                , PARAM_DOUBLE(1.)         },
  { "LLf"               , VAR(LLf)               , PARAM_DOUBLE(25.)        },
  { "LLz"               , VAR(LLz)               , PARAM_DOUBLE(400.*4)     },
  { "LLy"               , VAR(LLy)               , PARAM_DOUBLE(400.)       },

  { "background_n"      , VAR(background_n)      , PARAM_DOUBLE(.002)       },
  { "background_Te"     , VAR(background_Te)     , PARAM_DOUBLE(.001)       },
  { "background_Ti"     , VAR(background_Ti)     , PARAM_DOUBLE(.001)       },

  { "target_yl"         , VAR(target_yl)         , PARAM_DOUBLE(-100000.)   },
  { "target_yh"         , VAR(target_yh)         , PARAM_DOUBLE( 100000.)   },
  { "target_zwidth"     , VAR(target_zwidth)     , PARAM_DOUBLE(1.)         },

  { "no_initial_target" , VAR(no_initial_target) , PARAM_BOOL(false)        },

  { "heating_zl"        , VAR(heating_zl)        , PARAM_DOUBLE(-1.)        },
  { "heating_zh"        , VAR(heating_zh)        , PARAM_DOUBLE(1.)         },
  { "heating_xc"        , VAR(heating_xc)        , PARAM_DOUBLE(0.)         },
  { "heating_yc"        , VAR(heating_yc)        , PARAM_DOUBLE(0.)         },
  { "heating_rH"        , VAR(heating_rH)        , PARAM_DOUBLE(3.)         },

  { "LLs"               , VAR(LLs)               , MRC_VAR_DOUBLE           },
  { "LLn"               , VAR(LLn)               , MRC_VAR_DOUBLE           },
  { "target"            , VAR(target)            , MRC_VAR_OBJ(psc_target)  },
  { "inject"            , VAR(inject)            , MRC_VAR_OBJ(psc_inject)  },
  {},
};
#undef VAR

// ======================================================================
// PscHeatingSpotFoil

struct PscHeatingSpotFoilParams
{
  double zl; // in internal units (d_e)
  double zh;
  double xc;
  double yc;
  double rH;
  double T;
  double Mi;
};

struct PscHeatingSpotFoil : PscHeatingSpotFoilParams
{
  PscHeatingSpotFoil(const PscHeatingSpotFoilParams& params)
    : PscHeatingSpotFoilParams{params}
  {
    double width = zh - zl;
    fac = (8.f * pow(T, 1.5)) / (sqrt(Mi) * width);
    // FIXME, I don't understand the sqrt(Mi) in here
  }
  
  double operator()(const double *xx)
  {
    double x = xx[0], y = xx[1], z = xx[2];

    if (z <= zl || z >= zh) {
      return 0;
    }
    
    return fac * exp(-(sqr(x-xc) + sqr(y-yc)) / sqr(rH));
  }

private:
  double fac;
};

// ======================================================================
// PscFlatfoil
//
// eventually, a Psc replacement / derived class, but for now just
// pretending to be something like that
//
// things are missing from the generic step():
// - timing
// - psc_checks
// - pushp prep
// - marder

struct PscFlatfoil : Params
{
  using Mparticles_t = MparticlesDouble;
  using Mfields_t = MfieldsC;
#if 1 // generic_c
  using PushParticlesPusher_t = PushParticles__<Config2nd<dim_yz>>;
#else // 1vbec
  using PushParticlesPusher_t = PushParticles1vb<Config1vbec<Mparticles_t, Mfields_t, dim_yz>>;
#endif
  
  using Sort_t = SortCountsort2<Mparticles_t>;
  using Collision_t = Collision_<Mparticles_t, Mfields_t>;
  using PushFields_t = PushFields<Mfields_t>;
  using BndParticles_t = psc_bnd_particles_sub<Mparticles_t>;
  using Bnd_t = Bnd_<Mfields_t>;
  using BndFields_t = BndFieldsNone<Mfields_t>; // FIXME, why MfieldsC hardcoded???
  using Inject_t = Inject_<Mparticles_t, PscMfieldsC::sub_t>; // FIXME, shouldn't always use MfieldsC
  using Heating_t = Heating__<Mparticles_t>;
  using Balance_t = Balance_<PscMparticles<Mparticles_t>, PscMfields<Mfields_t>>;

  PscFlatfoil(psc *psc)
    : Params(psc->params),
      psc_{psc},
      sub_{psc_flatfoil(psc)},
      mprts_{dynamic_cast<Mparticles_t&>(*PscMparticlesBase{psc->particles}.sub())},
      mflds_{dynamic_cast<Mfields_t&>(*PscMfieldsBase{psc->flds}.sub())},
      collision_{psc_comm(psc), collision_interval, collision_nu},
      balance_{dynamic_cast<Balance_t&>(*PscBalanceBase{psc->balance}.sub())}
  {
    MPI_Comm comm = psc_comm(psc);
    
    bndp_.reset(new BndParticles_t{psc_->mrc_domain, psc_->grid()});
    bnd_.reset(new Bnd_t{psc_->grid(), psc_->mrc_domain, psc_->ibn});

    PscHeatingSpotFoilParams foil_params;
    foil_params.zl = sub_->heating_zl * sub_->d_i;
    foil_params.zh = sub_->heating_zh * sub_->d_i;
    foil_params.xc = sub_->heating_xc * sub_->d_i;
    foil_params.yc = sub_->heating_yc * sub_->d_i;
    foil_params.rH = sub_->heating_rH * sub_->d_i;
    foil_params.T  = .04;
    foil_params.Mi = sub_->heating_rH * psc->kinds[MY_ION].m;
    heating_.reset(new Heating_t{20, 0, 10000000, MY_ELECTRON, PscHeatingSpotFoil{foil_params}});

    bool inject_enable = true;
    bool inject_kind_n = MY_ELECTRON;
    bool inject_interval = 20;
    bool inject_tau = 40;
    inject_.reset(new Inject_t{comm, inject_enable, inject_interval, inject_tau, inject_kind_n,
	  sub_->target});
  }
  
  void step()
  {
    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = psc_comm(psc_);
    int timestep = psc_->timestep;

    balance_(psc_, mprts_);
    
    if (sort_interval > 0 && timestep % sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      sort_(mprts_);
    }
    
    if (collision_interval > 0 && ppsc->timestep % collision_interval == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      collision_(mprts_);
    }
    
    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    pushp_.push_mprts(mprts_, mflds_);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}
    
    // === field propagation B^{n+1/2} -> B^{n+1}
    pushf_.push_H<dim_yz>(mflds_, .5);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    (*bndp_)(mprts_);
    
    (*inject_)(mprts_);
    (*heating_)(mprts_);
    
    // === field propagation E^{n+1/2} -> E^{n+3/2}
    bndf_.fill_ghosts_H(mflds_);
    bnd_->fill_ghosts(mflds_, HX, HX + 3);
    
    bndf_.add_ghosts_J(mflds_);
    bnd_->add_ghosts(mflds_, JXI, JXI + 3);
    bnd_->fill_ghosts(mflds_, JXI, JXI + 3);
    
    pushf_.push_E<dim_yz>(mflds_, 1.);
    
    bndf_.fill_ghosts_E(mflds_);
    bnd_->fill_ghosts(mflds_, EX, EX + 3);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
    
    // === field propagation B^{n+1} -> B^{n+3/2}
    bndf_.fill_ghosts_E(mflds_);
    bnd_->fill_ghosts(mflds_, EX, EX + 3);
    
    pushf_.push_H<dim_yz>(mflds_, .5);
    
    bndf_.fill_ghosts_H(mflds_);
    bnd_->fill_ghosts(mflds_, HX, HX + 3);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}


    //psc_checks_continuity_after_particle_push(psc->checks, psc);

    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not it's natural time,
    // but div B should be == 0 at any time...)
    //psc_marder_run(psc->marder, psc->flds, psc->particles);
    
    //psc_checks_gauss(psc->checks, psc);

    //psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

  void setup()
  {
    setup_stats();
  }

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

      PscMparticlesBase mprts(psc_->particles);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
    
      psc_->timestep++; // FIXME, too hacky
      psc_output(psc_);

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts->get_n_prts();

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

  static void integrate(struct psc *psc)
  {
    PscFlatfoil flatfoil(psc);
    flatfoil.setup();
    flatfoil.integrate();
  }

private:
  psc* psc_;
  psc_flatfoil* sub_;
  Mparticles_t& mprts_;
  Mfields_t& mflds_;
  Balance_t& balance_;

  Sort_t sort_;
  Collision_t collision_;
  PushParticlesPusher_t pushp_;
  PushFields_t pushf_;
  std::unique_ptr<BndParticles_t> bndp_;
  std::unique_ptr<Bnd_t> bnd_;
  BndFields_t bndf_;

  std::unique_ptr<Heating_t> heating_;
  std::unique_ptr<Inject_t> inject_;
  
  int st_nr_particles;
  int st_time_step;
};

// ----------------------------------------------------------------------
// psc_flatfoil_create

static void
psc_flatfoil_create(struct psc *psc)
{
  struct psc_flatfoil *sub = psc_flatfoil(psc);

  psc_default_dimensionless(psc);

  psc->prm.nmax = 210001;
  psc->prm.nicell = 100;
  psc->prm.nr_populations = N_MY_KINDS;
  psc->prm.fractional_n_particles_per_cell = true;
  psc->prm.cfl = 0.75;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1600;
  psc->domain.gdims[2] = 1600*4;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");

  psc_target_set_type(sub->target, "foil");
}

// ----------------------------------------------------------------------
// psc_flatfoil_setup

static void
psc_flatfoil_setup(struct psc *psc)
{
  struct psc_flatfoil *sub = psc_flatfoil(psc);

  sub->LLs = 4. * sub->LLf;
  sub->LLn = .5 * sub->LLf;
  
  psc->domain.length[0] = 1.;
  psc->domain.length[1] = sub->LLy;
  psc->domain.length[2] = sub->LLz;

  // center around origin
  for (int d = 0; d < 3; d++) {
    psc->domain.corner[d] = -.5 * psc->domain.length[d];
  }

  // last population is neutralizing
  psc->kinds[MY_ELECTRON].q = -1.;
  psc->kinds[MY_ELECTRON].m = 1.;
  psc->kinds[MY_ELECTRON].name = strdup("e");

  psc->kinds[MY_ION     ].q = sub->Zi;
  psc->kinds[MY_ION     ].m = 100. * sub->Zi;  // FIXME, hardcoded mass ratio 100
  psc->kinds[MY_ION     ].name = strdup("i");

  sub->d_i = sqrt(psc->kinds[MY_ION].m / psc->kinds[MY_ION].q);

  psc_target_set_param_double(sub->target, "yl", sub->target_yl * sub->d_i);
  psc_target_set_param_double(sub->target, "yh", sub->target_yh * sub->d_i);
  psc_target_set_param_double(sub->target, "zl", - sub->target_zwidth * sub->d_i);
  psc_target_set_param_double(sub->target, "zh",   sub->target_zwidth * sub->d_i);

  psc_inject_set_param_int(sub->inject, "kind_n", MY_ELECTRON);
  psc_inject_set_param_obj(sub->inject, "target", sub->target);

  psc_setup_super(psc);
  psc_setup_member_objs_sub(psc);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "d_e = %g, d_i = %g\n", 1., sub->d_i);
  mpi_printf(comm, "lambda_De (background) = %g\n", sqrt(sub->background_Te));
}

// ----------------------------------------------------------------------
// psc_flatfoil_read

static void
psc_flatfoil_read(struct psc *psc, struct mrc_io *io)
{
  psc_read_super(psc, io);
}

// ----------------------------------------------------------------------
// psc_flatfoil_init_field

static double
psc_flatfoil_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_flatfoil *sub = psc_flatfoil(psc);

  double BB = sub->BB;

  switch (m) {
  case HY:
    return BB;

  default:
    return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_flatfoil_init_npt

static void
psc_flatfoil_init_npt(struct psc *psc, int pop, double x[3],
		      struct psc_particle_npt *npt)
{
  struct psc_flatfoil *sub = psc_flatfoil(psc);
  struct psc_target *target = sub->target;

  switch (pop) {
  case MY_ION:
    npt->n    = sub->background_n;
    npt->T[0] = sub->background_Ti;
    npt->T[1] = sub->background_Ti;
    npt->T[2] = sub->background_Ti;
    break;
  case MY_ELECTRON:
    npt->n    = sub->background_n;
    npt->T[0] = sub->background_Te;
    npt->T[1] = sub->background_Te;
    npt->T[2] = sub->background_Te;
    break;
  default:
    assert(0);
  }

  if (sub->no_initial_target) {
    return;
  }

  if (psc_target_is_inside(target, x)) {
    // replace values above by target values
    psc_target_init_npt(target, pop, x, npt);
  }
}

// ----------------------------------------------------------------------
// psc_ops "flatfoil"

struct psc_ops_flatfoil : psc_ops {
  psc_ops_flatfoil() {
    name             = "flatfoil";
    size             = sizeof(struct psc_flatfoil);
    param_descr      = psc_flatfoil_descr;
    create           = psc_flatfoil_create;
    setup            = psc_flatfoil_setup;
    read             = psc_flatfoil_read;
    init_field       = psc_flatfoil_init_field;
    init_npt         = psc_flatfoil_init_npt;
#ifdef USE_OWN_PSC_STEP
    integrate        = PscFlatfoil::integrate;
#endif
  }
} psc_flatfoil_ops;

// ======================================================================
// psc_target subclass "foil"

struct psc_target_foil {
  // params
  double yl;
  double yh;
  double zl;
  double zh;
  double n;
  double Te;
  double Ti;
};

#define psc_target_foil(target) mrc_to_subobj(target, struct psc_target_foil)

#define VAR(x) (void *)offsetof(struct psc_target_foil, x)
static struct param psc_target_foil_descr[] _mrc_unused = {
  { "yl"           , VAR(yl)           , PARAM_DOUBLE(0.)       },
  { "yh"           , VAR(yh)           , PARAM_DOUBLE(0.)       },
  { "zl"           , VAR(zl)           , PARAM_DOUBLE(0.)       },
  { "zh"           , VAR(zh)           , PARAM_DOUBLE(0.)       },
  { "n"            , VAR(n)            , PARAM_DOUBLE(1.)       },
  { "Te"           , VAR(Te)           , PARAM_DOUBLE(.001)     },
  { "Ti"           , VAR(Ti)           , PARAM_DOUBLE(.001)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_target_foil_is_inside

static bool
psc_target_foil_is_inside(struct psc_target *target, double x[3])
{
  struct psc_target_foil *sub = psc_target_foil(target);
  
  return (x[1] >= sub->yl && x[1] <= sub->yh &&
	  x[2] >= sub->zl && x[2] <= sub->zh);
}

// ----------------------------------------------------------------------
// psc_target_foil_init_npt

static void
psc_target_foil_init_npt(struct psc_target *target, int pop, double x[3],
			 struct psc_particle_npt *npt)
{
  struct psc_target_foil *sub = psc_target_foil(target);

  if (!psc_target_foil_is_inside(target, x)) {
    npt->n = 0;
    return;
  }

  switch (pop) {
  case MY_ION:
    npt->n    = sub->n;
    npt->T[0] = sub->Ti;
    npt->T[1] = sub->Ti;
    npt->T[2] = sub->Ti;
    break;
  case MY_ELECTRON:
    npt->n    = sub->n;
    npt->T[0] = sub->Te;
    npt->T[1] = sub->Te;
    npt->T[2] = sub->Te;
    break;
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_target "foil"

struct psc_target_ops_foil : psc_target_ops {
  psc_target_ops_foil() {
    name                = "foil";
    size                = sizeof(struct psc_target_foil);
    param_descr         = psc_target_foil_descr;
    is_inside           = psc_target_foil_is_inside;
    init_npt            = psc_target_foil_init_npt;
  }
} psc_target_ops_foil;

// ======================================================================
// main

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_target,
			      &psc_target_ops_foil);
  return psc_main(&argc, &argv, &psc_flatfoil_ops);
}

