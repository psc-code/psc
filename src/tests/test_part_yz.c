
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include "psc_bnd.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <mpi.h>

static bool do_dump = false;
static bool check_currents = true;

// ======================================================================
// psc_calc_rho

typedef double creal;

static inline creal
creal_sqrt(creal x)
{
  return sqrt(x);
}

static void
do_shift_particle_positions(particles_t *pp, double dt)
{
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;
    part->xi += vxi * dt;
    part->yi += vyi * dt;
    part->zi += vzi * dt;
  }
}

static void
psc_shift_particle_positions(struct psc *psc, mparticles_base_t *particles_base,
			     double dt)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  psc_foreach_patch(psc, p) {
    do_shift_particle_positions(psc_mparticles_get_patch(particles, p), dt);
  }

  psc_mparticles_put_cf(particles, particles_base);
}

// FIXME, this is pretty much the same as in calc_moments

static void
do_calc_rho_2nd(struct psc *psc, int p, particles_t *pp, fields_t *rho)
{
  creal fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  creal dxi = 1.f / psc->dx[0];
  creal dyi = 1.f / psc->dx[1];
  creal dzi = 1.f / psc->dx[2];

  struct psc_patch *patch = &psc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    creal xi = part->xi;
    creal yi = part->yi;
    creal zi = part->zi;
    creal u = (xi - patch->xb[0]) * dxi;
    creal v = (yi - patch->xb[1]) * dyi;
    creal w = (zi - patch->xb[2]) * dzi;
    int j1 = particle_real_nint(u);
    int j2 = particle_real_nint(v);
    int j3 = particle_real_nint(w);
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
		 mfields_base_t *rho_base)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);
  mfields_t *rho = psc_mfields_get_cf(rho_base, 0, 0);

  psc_mfields_zero(rho, 0);
  psc_foreach_patch(psc, p) {
    do_calc_rho_2nd(psc, p, psc_mparticles_get_patch(particles, p),
		    psc_mfields_get_patch(rho, p));
  }

  psc_mparticles_put_cf(particles, particles_base);
  psc_mfields_put_cf(rho, rho_base, 0, 1);

  psc_bnd_add_ghosts(psc->bnd, rho_base, 0, 1);
}

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

static struct psc_case *
create_test(const char *s_push_particles)
{
  struct psc_case *_case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, s_push_particles);
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_sort_set_param_int3(ppsc->sort, "blocksize", (int [3]) { 1, 8, 8 }); // FIXME
  psc_case_setup(_case);
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
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
  mfields_t *flds = psc_mfields_get_cf(flds_base, JXI, JXI + 3);
  mfields_t *div_j = psc_mfields_get_cf(div_j_base, 0, 0);

  psc_foreach_patch(psc, p) {
    do_calc_div_j(psc, p, psc_mfields_get_patch(flds, p), psc_mfields_get_patch(div_j, p));
  }

  psc_mfields_put_cf(flds, flds_base, 0, 0);
  psc_mfields_put_cf(div_j, div_j_base, 0, 1);

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

  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);
  psc_calc_rho_2nd(psc, particles, rho_m_base);
  psc_shift_particle_positions(psc, particles,  1. * psc->dt);
  psc_calc_rho_2nd(psc, particles, rho_p_base);
  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);

  // rho_p = (rho_p - rho_m) / dt
  psc_mfields_axpy(rho_p_base, -1., rho_m_base);
  psc_mfields_scale(rho_p_base, 1. / psc->dt);

  psc_calc_div_j(psc, flds, div_j_base);

  mfields_t *rho_p = psc_mfields_get_cf(rho_p_base, 0, 1);
  mfields_t *div_j = psc_mfields_get_cf(div_j_base, 0, 1);

  //  psc_dump_field(div_j, 0, "div_j");
  //  psc_dump_field(rho_p, 0, "dt_rho");

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

  psc_mfields_put_cf(rho_p, rho_p_base, 0, 0);
  psc_mfields_put_cf(div_j, div_j_base, 0, 0);

  //  psc_mfields_axpy(rho_p, +1., div_j);
  //  psc_dump_field(rho_p, 0, "cont_diff");

  assert(max_err <= eps);

  psc_mfields_destroy(rho_m_base);
  psc_mfields_destroy(rho_p_base);
  psc_mfields_destroy(div_j_base);
}

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields,
	 struct psc_case *(*create_test)(const char *), const char *push)
{
  printf("=== testing push_part_yz%s() %s %s\n", push, s_push_particles,
	 is_ref ? "(ref)" : "");

  struct psc_case *_case = create_test(s_push_particles);
  dump(s_push_particles, 0);
  if (strlen(push) == 0) {
    psc_push_particles_run(ppsc->push_particles, ppsc->particles, ppsc->flds);
  } else if (strcmp(push, "_a") == 0) {
    psc_push_particles_push_yz_a(ppsc->push_particles, ppsc->particles, ppsc->flds);
  } else if (strcmp(push, "_b") == 0) {
    psc_push_particles_push_yz_b(ppsc->push_particles, ppsc->particles, ppsc->flds);
  }
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  dump(s_push_particles, 1);
  if (strlen(push) == 0) { // only check continuity for full pusher
    psc_check_continuity(ppsc, ppsc->particles, ppsc->flds, 1e-14);
  }
  if (is_ref) {
    psc_save_particles_ref(ppsc, ppsc->particles);
    psc_save_fields_ref(ppsc, ppsc->flds);
  } else {
    psc_check_particles_ref(ppsc, ppsc->particles, eps_particles, "push_part_yz()");
    if (check_currents && strlen(push) == 0) { // only check currents for full pusher
      psc_check_currents_ref(ppsc, ppsc->flds, eps_fields);
    }
  }
  psc_case_destroy(_case);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_params_get_option_bool("dump", &do_dump);
  mrc_params_get_option_bool("check_currents", &check_currents);

  // ----------------------------------------------------------------------
  // push_yz_a

  run_test(true, "fortran", 0., 0., create_test, "_a");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "_a");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 1e-3, 1e-2, create_test, "_a");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_a");
#endif

  // ----------------------------------------------------------------------
  // push_yz_b

  run_test(true, "fortran", 0., 0., create_test, "_b");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "_b");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 2e-3, 1e-2, create_test, "_b");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_b");
#endif

  // ----------------------------------------------------------------------
  // push_yz

  run_test(true, "fortran", 0., 0., create_test, "");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 2e-3, 1e-3, create_test, "");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "");
#endif

  prof_print();

  MPI_Finalize();
}
