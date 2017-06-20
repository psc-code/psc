
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

struct psc_flatfoil {
  double BB;
  /* double TTe; */
  /* double TTi; */
  double n_bg;
  double Te_bg;
  double Ti_bg;
  double Zi;
  double LLf;
  double LLz;
  double LLy;

  // state
  double LLs;
  double LLn;
};

#define to_psc_flatfoil(psc) mrc_to_subobj(psc, struct psc_flatfoil)

#define VAR(x) (void *)offsetof(struct psc_flatfoil, x)
static struct param psc_flatfoil_descr[] = {
  { "BB"              , VAR(BB)              , PARAM_DOUBLE(.0)     },
  { "n_bg"            , VAR(n_bg)            , PARAM_DOUBLE(.002)   },
  { "Te_bg"           , VAR(Te_bg)           , PARAM_DOUBLE(.001)   },
  { "Ti_bg"           , VAR(Ti_bg)           , PARAM_DOUBLE(.001)   },
  { "Zi"              , VAR(Zi)              , PARAM_DOUBLE(1.)     },
  { "LLf"             , VAR(LLf)             , PARAM_DOUBLE(25.)    },
  { "LLz"             , VAR(LLz)             , PARAM_DOUBLE(400.*4) },
  { "LLy"             , VAR(LLy)             , PARAM_DOUBLE(400.)   },

  { "LLs"             , VAR(LLs)             , MRC_VAR_DOUBLE       },
  { "LLn"             , VAR(LLn)             , MRC_VAR_DOUBLE       },
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
}

// ----------------------------------------------------------------------
// psc_flatfoil_setup

static void
psc_flatfoil_setup(struct psc *psc)
{
  struct psc_flatfoil *sub = to_psc_flatfoil(psc);

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

  psc->kinds[MY_ION     ].q = sub->Zi;
  psc->kinds[MY_ION     ].m = 100. * sub->Zi;  // FIXME, hardcoded mass ratio 100

  psc_setup_super(psc);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "lambda_De_bg = %g\n", sqrt(sub->Te_bg));
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
  struct psc_flatfoil *sub = to_psc_flatfoil(psc);

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
  struct psc_flatfoil *sub = to_psc_flatfoil(psc);

  switch (pop) {
  case MY_ION:
    npt->n    = sub->n_bg;
    npt->T[0] = sub->Ti_bg;
    npt->T[1] = sub->Ti_bg;
    npt->T[2] = sub->Ti_bg;
    break;
  case MY_ELECTRON:
    npt->n    = sub->n_bg;
    npt->T[0] = sub->Te_bg;
    npt->T[1] = sub->Te_bg;
    npt->T[2] = sub->Te_bg;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_flatfoil_ops

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
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_flatfoil_ops);
}
