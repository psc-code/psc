
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_heating.h>
#include <psc_inject.h>
#include <psc_target_private.h>
#include <psc_event_generator_private.h>

#include <math.h>
#include <stdlib.h>

// ======================================================================
// psc subclass "flatfoil"

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
  struct psc_heating *heating;
  
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
  { "heating"           , VAR(heating)           , MRC_VAR_OBJ(psc_heating) },
  {},
};
#undef VAR

enum {
  MY_ION,
  MY_ELECTRON,
  N_MY_KINDS,
};

// ----------------------------------------------------------------------
// psc_flatfoil_create

static void
psc_flatfoil_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 210001;
  psc->prm.nicell = 100;
  psc->prm.gdims_in_terms_of_cells = true;
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

  psc_event_generator_set_type(psc->event_generator, "flatfoil");
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
  psc->kinds[MY_ELECTRON].name = "e";

  psc->kinds[MY_ION     ].q = sub->Zi;
  psc->kinds[MY_ION     ].m = 100. * sub->Zi;  // FIXME, hardcoded mass ratio 100
  psc->kinds[MY_ION     ].name = "i";

  sub->d_i = sqrt(psc->kinds[MY_ION].m / psc->kinds[MY_ION].q);

  psc_target_set_param_double(sub->target, "yl", sub->target_yl * sub->d_i);
  psc_target_set_param_double(sub->target, "yh", sub->target_yh * sub->d_i);
  psc_target_set_param_double(sub->target, "zl", - sub->target_zwidth * sub->d_i);
  psc_target_set_param_double(sub->target, "zh",   sub->target_zwidth * sub->d_i);
  psc_target_set_param_int(sub->target, "kind_electron", MY_ELECTRON);
  psc_target_set_param_int(sub->target, "kind_ion", MY_ION);

  psc_inject_set_param_int(sub->inject, "kind_n", MY_ELECTRON);
  psc_inject_set_param_obj(sub->inject, "target", sub->target);

  psc_heating_set_param_double(sub->heating, "zl", sub->heating_zl * sub->d_i);
  psc_heating_set_param_double(sub->heating, "zh", sub->heating_zh * sub->d_i);
  psc_heating_set_param_double(sub->heating, "xc", sub->heating_xc * sub->d_i);
  psc_heating_set_param_double(sub->heating, "yc", sub->heating_yc * sub->d_i);
  psc_heating_set_param_double(sub->heating, "rH", sub->heating_rH * sub->d_i);
  psc_heating_set_param_double(sub->heating, "Mi", psc->kinds[MY_ION].m);
  psc_heating_set_param_int(sub->heating, "kind", MY_ELECTRON);

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

struct psc_ops psc_flatfoil_ops = {
  .name             = "flatfoil",
  .size             = sizeof(struct psc_flatfoil),
  .param_descr      = psc_flatfoil_descr,
  .create           = psc_flatfoil_create,
  .setup            = psc_flatfoil_setup,
  .read             = psc_flatfoil_read,
  .init_field       = psc_flatfoil_init_field,
  .init_npt         = psc_flatfoil_init_npt,
};

// ======================================================================
// psc_event_generator subclass "flatfoil"

// ----------------------------------------------------------------------
// psc_event_generator_flatfoil_run

void
psc_event_generator_flatfoil_run(struct psc_event_generator *gen,
				 mparticles_base_t *mprts, mfields_base_t *mflds,
				 mphotons_t *mphotons)
{
  struct psc *psc = ppsc; // FIXME
  struct psc_flatfoil *sub = psc_flatfoil(psc);

  psc_inject_run(sub->inject, mprts, mflds);
  psc_heating_run(sub->heating, mprts, mflds);
}

// ----------------------------------------------------------------------
// psc_event_generator_ops "flatfoil"

struct psc_event_generator_ops psc_event_generator_flatfoil_ops = {
  .name                  = "flatfoil",
  .run                   = psc_event_generator_flatfoil_run,
};


// ======================================================================
// psc_target subclass "slab"

struct psc_target_slab {
  // params
  double yl;
  double yh;
  double zl;
  double zh;
  double n;
  double Te;
  double Ti;
  int kind_ion;
  int kind_electron;
};

#define psc_target_slab(target) mrc_to_subobj(target, struct psc_target_slab)

#define VAR(x) (void *)offsetof(struct psc_target_slab, x)
static struct param psc_target_slab_descr[] _mrc_unused = {
  { "yl"           , VAR(yl)           , PARAM_DOUBLE(0.)       },
  { "yh"           , VAR(yh)           , PARAM_DOUBLE(0.)       },
  { "zl"           , VAR(zl)           , PARAM_DOUBLE(0.)       },
  { "zh"           , VAR(zh)           , PARAM_DOUBLE(0.)       },
  { "n"            , VAR(n)            , PARAM_DOUBLE(1.)       },
  { "Te"           , VAR(Te)           , PARAM_DOUBLE(.001)     },
  { "Ti"           , VAR(Ti)           , PARAM_DOUBLE(.001)     },
  { "kind_ion"     , VAR(kind_ion)     , PARAM_INT(-1)          },
  { "kind_electron", VAR(kind_electron), PARAM_INT(-1)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_target_slab_is_inside

static bool
psc_target_slab_is_inside(struct psc_target *target, double x[3])
{
  struct psc_target_slab *sub = psc_target_slab(target);
  
  return (x[1] >= sub->yl && x[1] <= sub->yh &&
	  x[2] >= sub->zl && x[2] <= sub->zh);
}

// ----------------------------------------------------------------------
// psc_target_slab_init_npt

static void
psc_target_slab_init_npt(struct psc_target *target, int pop, double x[3],
			 struct psc_particle_npt *npt)
{
  struct psc_target_slab *sub = psc_target_slab(target);

  if (!psc_target_slab_is_inside(target, x)) {
    npt->n = 0;
    return;
  }

  if (pop == sub->kind_ion) {
    npt->n    = sub->n;
    npt->T[0] = sub->Ti;
    npt->T[1] = sub->Ti;
    npt->T[2] = sub->Ti;
  } else if (pop == sub->kind_electron) {
    npt->n    = sub->n;
    npt->T[0] = sub->Te;
    npt->T[1] = sub->Te;
    npt->T[2] = sub->Te;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_target "slab"

static struct psc_target_ops psc_target_ops_slab = {
  .name                = "slab",
  .size                = sizeof(struct psc_target_slab),
  .param_descr         = psc_target_slab_descr,
  .is_inside           = psc_target_slab_is_inside,
  .init_npt            = psc_target_slab_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_event_generator,
			      &psc_event_generator_flatfoil_ops);
  mrc_class_register_subclass(&mrc_class_psc_target,
			      &psc_target_ops_slab);
  return psc_main(&argc, &argv, &psc_flatfoil_ops);
}
