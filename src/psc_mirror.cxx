
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_diag.h>
#include <psc_diag_item_private.h>
#include <psc_fields_c.h>
#include "fields.hxx"

#include <mrc_params.h>

#include <math.h>

// ======================================================================

enum {
  MIRROR_ELECTRON,
  MIRROR_ION,
  MIRROR_HELIUM,
  NR_MIRROR_KINDS,
};

// ======================================================================

struct psc_mirror {
  double mi_over_me;
  double mh_over_mi;
  double vA_over_c;
  double beta_e_par;
  double beta_i_par;
  double beta_h_par;
  double Th_perp_over_Th_par;
  double Ti_perp_over_Ti_par;
  double Te_perp_over_Te_par;
  double theta_0;
  double B1;
  double heavy_ions_density;
};

#define to_psc_mirror(psc) mrc_to_subobj(psc, struct psc_mirror)

#define VAR(x) (void *)offsetof(struct psc_mirror, x)
static struct param psc_mirror_descr[] = {
  { "mi_over_me"         , VAR(mi_over_me)         , PARAM_DOUBLE(10.) },
  { "mh_over_mi"         , VAR(mh_over_mi)         , PARAM_DOUBLE(1.)  },
  { "vA_over_c"          , VAR(vA_over_c)          , PARAM_DOUBLE(.1) },
  { "beta_e_par"         , VAR(beta_e_par)         , PARAM_DOUBLE(.2)  },
  { "beta_i_par"         , VAR(beta_i_par)         , PARAM_DOUBLE(2.)  },
  { "beta_h_par"         , VAR(beta_h_par)         , PARAM_DOUBLE(0)  },
  { "Th_perp_over_Th_par", VAR(Th_perp_over_Th_par), PARAM_DOUBLE(0) },
  { "Ti_perp_over_Ti_par", VAR(Ti_perp_over_Ti_par), PARAM_DOUBLE(2.5) },
  { "Te_perp_over_Te_par", VAR(Te_perp_over_Te_par), PARAM_DOUBLE(1.)  },
  { "theta_0"            , VAR(theta_0)            , PARAM_DOUBLE(75.*M_PI/180) },
  { "B1"                 , VAR(B1)                 , PARAM_DOUBLE(0.)  },
  { "heavy_ions_density" , VAR(heavy_ions_density)  , PARAM_DOUBLE(0.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_diag_item_mirror

static void
psc_diag_item_mirror_run(struct psc_diag_item *item, struct psc *psc, double *result)
{
  struct psc_mirror *mirror = to_psc_mirror(psc);
  double B0 = mirror->vA_over_c;
  //  double HX0 = B0 * sin(mirror->theta_0);
  double HX0 = 0.;    
  double HY0 = 0.;
  double HZ0 = B0;
  //  double HZ0 = B0 * cos(mirror->theta_0);

  PscMfieldsC mf = psc->flds->get_as<PscMfieldsC>(HX, HX + 3);
  psc_foreach_patch(psc, p) {
    struct psc_patch *patch = &psc->patch[p];
    double fac = patch->dx[0] * patch->dx[1] * patch->dx[2];
    Fields3d<fields_c_t> F(mf[p]);
    psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
      result[0] += 
	(sqr(F(HX, ix,iy,iz) - HX0) +
	 sqr(F(HY, ix,iy,iz) - HY0) +
	 sqr(F(HZ, ix,iy,iz) - HZ0)) * fac;
    } foreach_3d_end;
  }
  mf.put_as(psc->flds, 0, 0);
}

struct psc_diag_item_ops psc_diag_item_mirror_ops = {
  .name = "mirror",
  .run = psc_diag_item_mirror_run,
  .nr_values = 1,
  .title = { "deltaH2" },
};

// ----------------------------------------------------------------------
// psc_mirror_create

static void
psc_mirror_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.nmax = 256000;
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;
  psc->prm.cfl = 0.98;

  psc->prm.nicell = 800;

  struct psc_kind kinds[NR_MIRROR_KINDS] = {
    [MIRROR_ELECTRON] = { .name = strdup("e"), .q = -1.,          },
    [MIRROR_ION]      = { .name = strdup("i"), .q =  1., .m = 1., },
    [MIRROR_HELIUM]   = { .name = strdup("h"), .q =  1.,          },
  };
  psc_set_kinds(psc, NR_MIRROR_KINDS, kinds);

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 512;

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

  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.; 
  psc->domain.length[2] = 52.; 

  struct psc_bnd_fields *bnd_fields = 
    psc_push_fields_get_bnd_fields(psc->push_fields);
  psc_bnd_fields_set_type(bnd_fields, "none");
}

// ----------------------------------------------------------------------
// psc_mirror_init_npt

static void
psc_mirror_init_npt(struct psc *psc, int kind, double x[3],
			    struct psc_particle_npt *npt)
{
  struct psc_mirror *mirror = to_psc_mirror(psc);

  double B0 = mirror->vA_over_c;
  double Te_par = mirror->beta_e_par * sqr(B0) / 2.;
  double Te_perp = mirror->Te_perp_over_Te_par * Te_par;
  double Ti_par = mirror->beta_i_par * sqr(B0) / 2.;
  double Ti_perp = mirror->Ti_perp_over_Ti_par * Ti_par;
  double Th_par = mirror->beta_h_par * sqr(B0) / 2.;
  double Th_perp = mirror->Th_perp_over_Th_par * Th_par;
  double n_h = mirror->heavy_ions_density;

  switch (kind) {
  case MIRROR_ELECTRON: // electrons
    npt->n = 1.;
    npt->q = -1.;
    // npt->m = 1. / mirror->mi_over_me;

    npt->T[0] = Te_perp;
    npt->T[1] = Te_perp;
    npt->T[2] = Te_par;
    break;
  case MIRROR_ION: // ions
    npt->n = 1. - n_h;
    npt->q = 1.;
    // npt->m = 1.;

    npt->T[0] = Ti_perp;
    npt->T[1] = Ti_perp;
    npt->T[2] = Ti_par;
    break;
  case MIRROR_HELIUM: // helium
    npt->n = n_h;
    npt->q = 1.;
    // npt->m = mirror->mh_over_mi;

    npt->T[0] = Th_perp;
    npt->T[1] = Th_perp;
    npt->T[2] = Th_par;
    break;
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_mirror_init_field

static double
psc_mirror_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_mirror *mirror = to_psc_mirror(psc);
  double B0 = mirror->vA_over_c;
  double B1 = mirror->B1; 
  double kz = (32. * M_PI / psc->domain.length[2]) * cos(mirror->theta_0);
  double ky = (8. * M_PI / psc->domain.length[1]) * sin(mirror->theta_0);
  double z = x[2];
  double y = x[1];

  switch (m) {
  case HZ:
    return B0 + B1 * cos(ky*y + kz*z);

  default:
    return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_mirror_setup

static void
psc_mirror_setup(struct psc *psc)
{
  struct psc_mirror *mirror = to_psc_mirror(psc);

  double B0 = mirror->vA_over_c;
  double Te_par = mirror->beta_e_par * sqr(B0) / 2.;
  double Te_perp = mirror->Te_perp_over_Te_par * Te_par; 
  double Ti_par = mirror->beta_i_par * sqr(B0) / 2.;
  double Ti_perp = mirror->Ti_perp_over_Ti_par * Ti_par;
  double Th_par = mirror->beta_h_par * sqr(B0) / 2.;
  double Th_perp = mirror->Th_perp_over_Th_par * Th_par;
  double me = 1. / mirror->mi_over_me;
  double mh = mirror->mh_over_mi;

  // set particle kind parameters
  assert(psc->nr_kinds == NR_MIRROR_KINDS);
  psc->kinds[MIRROR_ELECTRON].m = me;
  psc->kinds[MIRROR_ION     ].m = 1.;
  psc->kinds[MIRROR_HELIUM  ].m = mh;

  psc_setup_super(psc);

  MPI_Comm comm = psc_comm(psc);
  mpi_printf(comm, "d_i = 1., d_e = %g, d_h = %g\n", sqrt(me), sqrt(mh));
  mpi_printf(comm, "om_ci = %g, om_ce = %g, om_ch = %g\n", B0, B0 / me, B0 / mh);
  mpi_printf(comm, "\n");
  mpi_printf(comm, "v_i,perp = %g [c] T_i,perp = %g\n", sqrt(2*Ti_perp), Ti_perp);
  mpi_printf(comm, "v_i,par  = %g [c] T_i,par = %g\n", sqrt(Ti_par), Ti_par);
  mpi_printf(comm, "v_e,perp = %g [c] T_e,perp = %g\n", sqrt(2*Te_perp / me), Te_perp);
  mpi_printf(comm, "v_e,par  = %g [c] T_e,par = %g\n", sqrt(Te_par / me), Te_par);
  mpi_printf(comm, "v_h,perp = %g [c] T_h,perp = %g\n", sqrt(2*Th_perp / mh), Th_perp);
  mpi_printf(comm, "v_h,par  = %g [c] T_h,par = %g\n", sqrt(Th_par / mh), Th_par);
  mpi_printf(comm, "\n");
}

// ======================================================================
// psc_mirror_ops

struct psc_ops psc_mirror_ops = {
  .name             = "mirror",
  .size             = sizeof(struct psc_mirror),
  .param_descr      = psc_mirror_descr,
  .create           = psc_mirror_create,
  .setup            = psc_mirror_setup,
  .init_npt         = psc_mirror_init_npt,
  .init_field       = psc_mirror_init_field
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_diag_item,
                              &psc_diag_item_mirror_ops);
  return psc_main(&argc, &argv, &psc_mirror_ops);
}

