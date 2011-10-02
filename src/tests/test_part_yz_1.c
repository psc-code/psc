
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_sort.h"
#include "psc_randomize.h"
#include "psc_bnd.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

static unsigned int mask;

static bool do_dump = false;
static bool check_currents = true;
static bool check_particles = true;

// ======================================================================
// psc_calc_rho

typedef double creal;

static inline creal
creal_sqrt(creal x)
{
  return sqrt(x);
}

// FIXME, this is pretty much the same as in calc_moments
// FIXME, check continuity in _yz() (and others), too

static void
do_calc_rho_2nd(struct psc *psc, int p, particles_t *pp,
		fields_t *rho, double dt)
{
  creal fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  creal dxi = 1.f / psc->dx[0];
  creal dyi = 1.f / psc->dx[1];
  creal dzi = 1.f / psc->dx[2];

  struct psc_patch *patch = &psc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;
    creal xi = part->xi + vxi * dt;
    creal yi = part->yi + vyi * dt;
    creal zi = part->zi + vzi * dt;
    creal u = (xi - patch->xb[0]) * dxi;
    creal v = (yi - patch->xb[1]) * dyi;
    creal w = (zi - patch->xb[2]) * dzi;
    int j1 = particle_base_real_nint(u);
    int j2 = particle_base_real_nint(v);
    int j3 = particle_base_real_nint(w);
    creal h1 = j1-u;
    creal h2 = j2-v;
    creal h3 = j3-w;
      
    creal gmx=.5f*(.5f+h1)*(.5f+h1);
    creal gmy=.5f*(.5f+h2)*(.5f+h2);
    creal gmz=.5f*(.5f+h3)*(.5f+h3);
    creal g0x=.75f-h1*h1;
    creal g0y=.75f-h2*h2;
    creal g0z=.75f-h3*h3;
    creal g1x=.5f*(.5f-h1)*(.5f-h1);
    creal g1y=.5f*(.5f-h2)*(.5f-h2);
    creal g1z=.5f*(.5f-h3)*(.5f-h3);
      
    if (psc->domain.gdims[0] == 1) {
      j1 = 0; gmx = 0.; g0x = 1.; g1x = 0.;
    }
    if (psc->domain.gdims[1] == 1) {
      j2 = 0; gmy = 0.; g0y = 1.; g1y = 0.;
    }
    if (psc->domain.gdims[2] == 1) {
      j3 = 0; gmz = 0.; g0z = 1.; g1z = 0.;
    }

    creal fnq = part->qni * part->wni * fnqs;
    F3(rho,0, j1-1,j2-1,j3-1) += fnq*gmx*gmy*gmz;
    F3(rho,0, j1  ,j2-1,j3-1) += fnq*g0x*gmy*gmz;
    F3(rho,0, j1+1,j2-1,j3-1) += fnq*g1x*gmy*gmz;
    F3(rho,0, j1-1,j2  ,j3-1) += fnq*gmx*g0y*gmz;
    F3(rho,0, j1  ,j2  ,j3-1) += fnq*g0x*g0y*gmz;
    F3(rho,0, j1+1,j2  ,j3-1) += fnq*g1x*g0y*gmz;
    F3(rho,0, j1-1,j2+1,j3-1) += fnq*gmx*g1y*gmz;
    F3(rho,0, j1  ,j2+1,j3-1) += fnq*g0x*g1y*gmz;
    F3(rho,0, j1+1,j2+1,j3-1) += fnq*g1x*g1y*gmz;
    F3(rho,0, j1-1,j2-1,j3  ) += fnq*gmx*gmy*g0z;
    F3(rho,0, j1  ,j2-1,j3  ) += fnq*g0x*gmy*g0z;
    F3(rho,0, j1+1,j2-1,j3  ) += fnq*g1x*gmy*g0z;
    F3(rho,0, j1-1,j2  ,j3  ) += fnq*gmx*g0y*g0z;
    F3(rho,0, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3(rho,0, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3(rho,0, j1-1,j2+1,j3  ) += fnq*gmx*g1y*g0z;
    F3(rho,0, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3(rho,0, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3(rho,0, j1-1,j2-1,j3+1) += fnq*gmx*gmy*g1z;
    F3(rho,0, j1  ,j2-1,j3+1) += fnq*g0x*gmy*g1z;
    F3(rho,0, j1+1,j2-1,j3+1) += fnq*g1x*gmy*g1z;
    F3(rho,0, j1-1,j2  ,j3+1) += fnq*gmx*g0y*g1z;
    F3(rho,0, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3(rho,0, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3(rho,0, j1-1,j2+1,j3+1) += fnq*gmx*g1y*g1z;
    F3(rho,0, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3(rho,0, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

void
psc_calc_rho_2nd(struct psc *psc, mparticles_base_t *particles_base,
		 mfields_base_t *rho_base, double dt)
{
  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);
  mfields_t *rho = psc_mfields_get_from(0, 0, rho_base);

  psc_mfields_zero(rho, 0);
  psc_foreach_patch(psc, p) {
    do_calc_rho_2nd(psc, p, &particles.p[p], psc_mfields_get_patch(rho, p), dt);
  }

  psc_mparticles_put_to(&particles, particles_base);
  psc_mfields_put_to(rho, 0, 1, rho_base);

  psc_bnd_add_ghosts(psc->bnd, rho_base, 0, 1);
}

// ======================================================================

static void
do_calc_rho_1st(struct psc *psc, int p, particles_t *pp,
		fields_t *rho, double dt)
{
  creal fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  creal dxi = 1.f / psc->dx[0];
  creal dyi = 1.f / psc->dx[1];
  creal dzi = 1.f / psc->dx[2];

  struct psc_patch *patch = &psc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;
    creal xi = part->xi + vxi * dt;
    creal yi = part->yi + vyi * dt;
    creal zi = part->zi + vzi * dt;
    creal u = (xi - patch->xb[0]) * dxi;
    creal v = (yi - patch->xb[1]) * dyi;
    creal w = (zi - patch->xb[2]) * dzi;
    int j1 = particle_base_real_fint(u);
    int j2 = particle_base_real_fint(v);
    int j3 = particle_base_real_fint(w);
    creal h1 = u-j1;
    creal h2 = v-j2;
    creal h3 = w-j3;
      
    creal g0x=1.f - h1;
    creal g0y=1.f - h2;
    creal g0z=1.f - h3;
    creal g1x=h1;
    creal g1y=h2;
    creal g1z=h3;
      
    if (psc->domain.gdims[0] == 1) {
      j1 = 0; g0x = 1.; g1x = 0.;
    }
    if (psc->domain.gdims[1] == 1) {
      j2 = 0; g0y = 1.; g1y = 0.;
    }
    if (psc->domain.gdims[2] == 1) {
      j3 = 0; g0z = 1.; g1z = 0.;
    }

    creal fnq = part->qni * part->wni * fnqs;
    F3(rho,0, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3(rho,0, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3(rho,0, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3(rho,0, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3(rho,0, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3(rho,0, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3(rho,0, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3(rho,0, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

void
psc_calc_rho_1st(struct psc *psc, mparticles_base_t *particles_base,
		 mfields_base_t *rho_base, double dt)
{
  mparticles_t particles;
  psc_mparticles_get_from(&particles, particles_base);
  mfields_t *rho = psc_mfields_get_from(0, 0, rho_base);

  psc_mfields_zero(rho, 0);
  psc_foreach_patch(psc, p) {
    do_calc_rho_1st(psc, p, &particles.p[p], psc_mfields_get_patch(rho, p), dt);
  }

  psc_mparticles_put_to(&particles, particles_base);
  psc_mfields_put_to(rho, 0, 1, rho_base);

  psc_bnd_add_ghosts(psc->bnd, rho_base, 0, 1);
}

// ======================================================================
// psc_calc_div_j

static void
do_calc_div_j(struct psc *psc, int p, fields_t *flds, fields_t *div_j)
{
  creal h[3];
  for (int d = 0; d < 3; d++) {
    if (psc->domain.gdims[d] == 1) {
      h[d] = 0.;
    } else {
      h[d] = 1. / psc->dx[d];
    }
  }

  psc_foreach_3d_g(psc, p, jx, jy, jz) {
    F3(div_j,0, jx,jy,jz) =
      (F3(flds,JXI, jx,jy,jz) - F3(flds,JXI, jx-1,jy,jz)) * h[0] +
      (F3(flds,JYI, jx,jy,jz) - F3(flds,JYI, jx,jy-1,jz)) * h[1] +
      (F3(flds,JZI, jx,jy,jz) - F3(flds,JZI, jx,jy,jz-1)) * h[2];
  } psc_foreach_3d_g_end;
}

static void
psc_calc_div_j(struct psc *psc, mfields_base_t *flds_base, mfields_base_t *div_j_base)
{
  mfields_t *flds = psc_mfields_get_from(JXI, JXI + 3, flds_base);
  mfields_t *div_j = psc_mfields_get_from(0, 0, div_j_base);

  psc_foreach_patch(psc, p) {
    do_calc_div_j(psc, p, psc_mfields_get_patch(flds, p), psc_mfields_get_patch(div_j, p));
  }

  psc_mfields_put_to(flds, 0, 0, flds_base);
  psc_mfields_put_to(div_j, 0, 1, div_j_base);

  psc_bnd_add_ghosts(psc->bnd, div_j_base, 0, 1);
}

// ======================================================================
// psc_check_continuity

static mfields_base_t *
fld_create(struct psc *psc)
{
  mfields_base_t *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_setup(fld);

  return fld;
}

static void
psc_check_continuity(struct psc *psc, mparticles_base_t *particles,
		     mfields_base_t *flds, double eps)
{
  mfields_base_t *rho_m_base = fld_create(psc);
  mfields_base_t *rho_p_base = fld_create(psc);
  mfields_base_t *div_j_base = fld_create(psc);

  psc_calc_rho_1st(psc, particles, rho_m_base, -.5 * psc->dt);
  psc_calc_rho_1st(psc, particles, rho_p_base,  .5 * psc->dt);

  psc_mfields_axpy(rho_p_base, -1., rho_m_base);
  psc_mfields_scale(rho_p_base, 1. / psc->dt);

  psc_calc_div_j(psc, flds, div_j_base);

  //  psc_dump_field(div_j, 0, "div_j");
  //  psc_dump_field(rho_p, 0, "dt_rho");

  mfields_t *rho_p = psc_mfields_get_from(0, 1, rho_p_base);
  mfields_t *div_j = psc_mfields_get_from(0, 1, div_j_base);

  double max_err = 0.;
  psc_foreach_patch(psc, p) {
    fields_t *p_rho_p = psc_mfields_get_patch(rho_p, p);
    fields_t *p_div_j = psc_mfields_get_patch(div_j, p);
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      creal dt_rho = F3(p_rho_p,0, jx,jy,jz);
      creal div_j = F3(p_div_j,0, jx,jy,jz);
      max_err = fmax(max_err, fabs(dt_rho + div_j));
      if (fabs(dt_rho + div_j) > eps) {
	printf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
	       dt_rho, -div_j, dt_rho + div_j);
      }
    } psc_foreach_3d_end;
  }
  printf("continuity: max_err = %g (thres %g)\n", max_err, eps);

  psc_mfields_put_to(rho_p, 0, 0, rho_p_base);
  psc_mfields_put_to(div_j, 0, 0, div_j_base);

  //  psc_mfields_base_axpy(rho_p, +1., div_j);
  //  psc_dump_field(rho_p, 0, "cont_diff");

  assert(max_err <= eps);

  psc_mfields_destroy(rho_m_base);
  psc_mfields_destroy(rho_p_base);
  psc_mfields_destroy(div_j_base);
}


// ======================================================================
// psc_test

struct psc_test {
};

#define to_psc_test(psc) mrc_to_subobj(psc, struct psc_test)

// ----------------------------------------------------------------------
// psc_test_create

static void
psc_test_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 100;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 16;
  psc->domain.gdims[2] = 16;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  psc_sort_set_type(psc->sort, "countsort2");
#ifdef USE_CUDA // FIXME
#if BLOCKSIZE_X == 1 && BLOCKSIZE_Y == 4 && BLOCKSIZE_Z == 4
  psc_sort_set_param_int3(ppsc->sort, "blocksize", (int [3]) { 1, 8, 8 });
#else
#error TBD
#endif
#endif
}

// ----------------------------------------------------------------------
// psc_test_set_from_options

static void
psc_test_set_from_options(struct psc *psc)
{
  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;
}

// ----------------------------------------------------------------------
// psc_test_init_field

static double
psc_test_init_field(struct psc *psc, double x[3], int m)
{
#if 0
  switch (m) {
  case EY: return 1.;
  default: return 0.;
  }
#endif

  switch (m) {
  case EX: return x[1] + x[2];
  case EY: return x[1] + x[2];
  case EZ: return x[1] + x[2];
  case HX: return x[1] + x[2];
  case HY: return x[1] + x[2];
  case HZ: return x[1] + x[2];
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_init_npt

static void
psc_test_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  npt->n = 1.;
  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    break;
  default:
    assert(0);
  }
}

static void
psc_test_step(struct psc *psc)
{
  psc_output(psc);

  // field propagation n*dt -> (n+0.5)*dt
  psc_push_fields_step_a(psc->push_fields, psc->flds);

  // particle propagation n*dt -> (n+1.0)*dt
  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);

  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);

  // field propagation (n+0.5)*dt -> (n+1.0)*dt
  psc_push_fields_step_b(psc->push_fields, psc->flds);
}

// ======================================================================
// psc_test_ops

struct psc_ops psc_test_ops = {
  .name             = "test",
  .size             = sizeof(struct psc_test),
  .create           = psc_test_create,
  .set_from_options = psc_test_set_from_options,
  .init_field       = psc_test_init_field,
  .init_npt         = psc_test_init_npt,
  .step             = psc_test_step,
};

static void
dump(const char *basename, int cnt)
{
  if (!do_dump)
    return;

  char s[200];
  sprintf(s, "part_%s_%d", basename, cnt);
  psc_dump_particles(ppsc->particles, s);
  sprintf(s, "jx_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JXI, s);
  sprintf(s, "jy_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JYI, s);
  sprintf(s, "jz_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JZI, s);
}

static struct psc *
create_test(const char *s_push_particles)
{
  struct psc *psc = psc_create(MPI_COMM_WORLD);
  psc_set_from_options(psc);
  psc_push_particles_set_type(psc->push_particles, s_push_particles);
  psc_sort_set_type(psc->sort, "countsort2");
  psc_randomize_set_type(psc->randomize, "c");
  psc_setup(psc);
#if 0
  psc->particles->p[0].particles[0] = psc->particles->p[0].particles[1];
  psc->particles->p[0].particles[0].yi = 2./16;
  psc->particles->p[0].particles[0].zi = 3.5/16;
  psc->particles->p[0].n_part = 1;
#endif
  psc_sort_set_param_int(ppsc->sort, "mask", mask);
  psc_randomize_run(psc->randomize, psc->particles);
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  psc_sort_run(psc->sort, psc->particles);

  return psc;
}

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields,
	 struct psc *(*create_test)(const char *), const char *push)
{
  printf("=== testing push_part_yz%s() %s %s\n", push, s_push_particles,
	 is_ref ? "(ref)" : "");

  struct psc *psc = create_test(s_push_particles);
  dump(s_push_particles, 0);
  if (strlen(push) == 0) {
    psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  } else if (strcmp(push, "_a") == 0) {
    psc_push_particles_push_yz_a(psc->push_particles, psc->particles, psc->flds);
  } else if (strcmp(push, "_b") == 0) {
    psc_push_particles_push_yz_b(psc->push_particles, psc->particles, psc->flds);
  }
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  psc_sort_run(psc->sort, psc->particles);
  dump(s_push_particles, 1);
  if (is_ref) {
    psc_save_particles_ref(psc, psc->particles);
    psc_save_fields_ref(psc, psc->flds);
  } else {
    if (check_particles) {
      psc_check_particles_ref(psc, psc->particles, eps_particles, "push_part_yz()");
    }
    if (check_currents && strlen(push) == 0) { // only check currents for full pusher
      psc_check_currents_ref(psc, psc->flds, eps_fields);
      psc_check_continuity(psc, psc->particles, psc->flds, eps_fields);
    }
  }
  psc_destroy(psc);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  mrc_params_get_option_bool("dump", &do_dump);
  mrc_params_get_option_bool("check_currents", &check_currents);
  mrc_params_get_option_bool("check_particles", &check_particles);

  // ----------------------------------------------------------------------
  // push_yz 1st order

  run_test(true, "1st", 0., 0., create_test, "");
  // run again to check continuity
  run_test(false, "1st", 1e-7, 1e-7, create_test, "");

  // since the fields are linear functions of position, 1st order / 2nd order
  // field interpolation should give the same result
  run_test(false, "generic_c", 1e-7, 1e-0, create_test, "");

#ifdef USE_CUDA
  run_test(false, "cuda_1st", 1e-6, 1e-3, create_test, "");
#endif

  // ----------------------------------------------------------------------
  // push_yz 1st order vb

  mask = 15;

  run_test(true, "1vb", 0., 0., create_test, "");
  // run again to check continuity
  run_test(false, "1vb", 1e-7, 1e-3, create_test, "");

#ifdef USE_CUDA
  run_test(false, "cuda_1vb", 1e-6, 1e-3, create_test, "");
#endif

  prof_print();

  MPI_Finalize();
}
