

#include "psc.h"
#include "psc_glue.h"

#include "psc_push_fields.h"
#include "psc_bnd_fields_private.h"
#include "psc_pulse.h"
#include "psc_fields_c.h"

#include <assert.h>

#define PSC_set_globals_F77 F77_FUNC_(psc_set_globals, PSC_SET_GLOBALS)
#define PSC_set_patch_F77 F77_FUNC_(psc_set_patch, PSC_SET_PATCH)
#define PSC_set_params_F77 F77_FUNC_(psc_set_params, PSC_SET_PARAMS)
#define PSC_set_domain_F77 F77_FUNC_(psc_set_domain, PSC_SET_DOMAIN)
#define PSC_set_coeff_F77 F77_FUNC_(psc_set_coeff, PSC_SET_COEFF)
#define PSC_set_timestep_F77 F77_FUNC(psc_set_timestep, PSC_SET_TIMESTEP)
#define OUT_params_set_F77 F77_FUNC_(out_params_set, OUT_PARAMS_SET)
#define SETUP_field_F77 F77_FUNC_(setup_field, SETUP_FIELD)
#define PIC_push_part_xy_F77 F77_FUNC(pic_push_part_xy,PIC_PUSH_PART_XY)
#define PIC_push_part_xz_F77 F77_FUNC(pic_push_part_xz,PIC_PUSH_PART_XZ)
#define PIC_push_part_yz_F77 F77_FUNC(pic_push_part_yz,PIC_PUSH_PART_YZ)
#define PIC_push_part_xyz_F77 F77_FUNC(pic_push_part_xyz,PIC_PUSH_PART_XYZ)
#define PIC_push_part_z_F77 F77_FUNC(pic_push_part_z,PIC_PUSH_PART_Z)
#define PIC_push_part_z_vay_F77 F77_FUNC(pic_push_part_z_vay,PIC_PUSH_PART_Z_VAY)
#define PIC_randomize_F77 F77_FUNC(pic_randomize,PIC_RANDOMIZE)
#define PIC_bin_coll_F77 F77_FUNC(pic_bin_coll,PIC_BIN_COLL)
#define PIC_find_cell_indices_F77 F77_FUNC(pic_find_cell_indices,PIC_FIND_CELL_INDICES)
#define SET_param_pml_F77 F77_FUNC(set_param_pml,SET_PARAM_PML)
#define PIC_msa_e_F77 F77_FUNC_(pic_msa_e, PIC_MSA_E)
#define PIC_msa_h_F77 F77_FUNC_(pic_msa_h, PIC_MSA_H)
#define PIC_msb_h_F77 F77_FUNC_(pic_msb_h, PIC_MSB_H)
#define PIC_msb_e_F77 F77_FUNC_(pic_msb_e, PIC_MSB_E)
#define PIC_pml_msa_F77 F77_FUNC_(pic_pml_msa, PIC_PML_MSA)
#define PIC_pml_msb_F77 F77_FUNC_(pic_pml_msb, PIC_PML_MSB)
#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)
#define PIC_fill_ghosts_h_b_F77 F77_FUNC_(pic_fill_ghosts_h_b, PIC_FILL_GHOSTS_H_B)
#define PML_coeff_init_param_F77 F77_FUNC_(pml_coeff_init_param, PML_COEFF_INIT_PARAM)

#define C_p_pulse_x1_F77 F77_FUNC(c_p_pulse_x1,C_P_PULSE_X1)
#define C_s_pulse_x1_F77 F77_FUNC(c_s_pulse_x1,C_S_PULSE_X1)
#define C_p_pulse_x2_F77 F77_FUNC(c_p_pulse_x2,C_P_PULSE_X2)
#define C_s_pulse_x2_F77 F77_FUNC(c_s_pulse_x2,C_S_PULSE_X2)

#define C_p_pulse_y1_F77 F77_FUNC(c_p_pulse_y1,C_P_PULSE_Y1)
#define C_s_pulse_y1_F77 F77_FUNC(c_s_pulse_y1,C_S_PULSE_Y1)
#define C_p_pulse_y2_F77 F77_FUNC(c_p_pulse_y2,C_P_PULSE_Y2)
#define C_s_pulse_y2_F77 F77_FUNC(c_s_pulse_y2,C_S_PULSE_Y2)

#define C_p_pulse_z1_F77 F77_FUNC(c_p_pulse_z1,C_P_PULSE_Z1)
#define C_s_pulse_z1_F77 F77_FUNC(c_s_pulse_z1,C_S_PULSE_Z1)
#define C_p_pulse_z2_F77 F77_FUNC(c_p_pulse_z2,C_P_PULSE_Z2)
#define C_s_pulse_z2_F77 F77_FUNC(c_s_pulse_z2,C_S_PULSE_Z2)

void PSC_set_globals_F77(f_real *cori, f_real *alpha, f_real *eta);

void PSC_set_patch_F77(f_int *imn, f_int *imx, f_int *rd,
		       f_real *dt, f_real *dxyz);

void PSC_set_params_F77(f_real *qq, f_real *mm, f_real *tt, f_real *cc, f_real *eps0,
			f_int *nmax, f_real *lw, f_real *i0, f_real *n0,
			f_real *e0, f_real *b0, f_real *j0, f_real *rho0, f_real *phi0,
			f_real *a0);

void PSC_set_domain_F77(f_real *length, f_int *itot, f_int *in, f_int *ix,
			f_int *bnd_fld_lo, f_int *bnd_fld_hi, f_int *bnd_part,
			f_int *nproc, f_int *nghost, f_int *use_pml);

void PSC_set_coeff_F77(f_real *beta,
		       f_real *wl, f_real *ld, f_real *vos, f_real *vt, f_real *wp);

void PSC_set_timestep_F77(f_int *n);

void SETUP_field_F77(void);

void OUT_params_set_F77(f_int *np, f_int *nnp);

void PIC_push_part_xy_F77(f_int *niloc, particle_fortran_t *p_niloc,
			  f_real *p2A, f_real *p2B,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_xz_F77(f_int *niloc, particle_fortran_t *p_niloc,
			  f_real *p2A, f_real *p2B,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_yz_F77(f_int *niloc, particle_fortran_t *p_niloc,
			  f_real *p2A, f_real *p2B,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_xyz_F77(f_int *niloc, particle_fortran_t *p_niloc,
		 	   f_real *p2A, f_real *p2B,
			   f_real *jxi, f_real *jyi, f_real *jzi,
			   f_real *ex, f_real *ey, f_real *ez,
			   f_real *bx, f_real *by, f_real *bz);

void PIC_push_part_z_F77(f_int *niloc, particle_fortran_t *p_niloc,
			 f_real *p2A, f_real *p2B,
			 f_real *jxi, f_real *jyi, f_real *jzi,
			 f_real *ex, f_real *ey, f_real *ez,
			 f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_z_vay_F77(f_int *niloc, particle_fortran_t *p_niloc,
			     f_real *p2A, f_real *p2B,
			     f_real *jxi, f_real *jyi, f_real *jzi,
			     f_real *ex, f_real *ey, f_real *ez,
			     f_real *hx, f_real *hy, f_real *hz);

void PIC_randomize_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_bin_coll_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_find_cell_indices_F77(f_int *niloc, particle_fortran_t *p_niloc);
void SET_param_pml_F77(f_int *thick, f_int *cushion, f_int *size, f_int *order);

void PIC_msa_e_F77(f_real *ex, f_real *ey, f_real *ez,
		   f_real *hx, f_real *hy, f_real *hz,
		   f_real *jxi, f_real *jyi, f_real *jzi);
void PIC_msa_h_F77(f_real *ex, f_real *ey, f_real *ez,
		   f_real *hx, f_real *hy, f_real *hz,
		   f_real *jxi, f_real *jyi, f_real *jzi);
void PIC_msb_h_F77(f_real *ex, f_real *ey, f_real *ez,
		   f_real *hx, f_real *hy, f_real *hz,
		   f_real *jxi, f_real *jyi, f_real *jzi);
void PIC_msb_e_F77(f_real *ex, f_real *ey, f_real *ez,
		   f_real *hx, f_real *hy, f_real *hz,
		   f_real *jxi, f_real *jyi, f_real *jzi);
void PIC_pml_msa_F77(f_real *ex, f_real *ey, f_real *ez,
		     f_real *hx, f_real *hy, f_real *hz,
		     f_real *dx, f_real *dy, f_real *dz,
		     f_real *bx, f_real *by, f_real *bz,
		     f_real *jxi, f_real *jyi, f_real *jzi,
		     f_real *eps, f_real *mu);
void PIC_pml_msb_F77(f_real *ex, f_real *ey, f_real *ez,
		     f_real *hx, f_real *hy, f_real *hz,
		     f_real *dx, f_real *dy, f_real *dz,
		     f_real *bx, f_real *by, f_real *bz,
		     f_real *jxi, f_real *jyi, f_real *jzi,
		     f_real *eps, f_real *mu);
void PIC_fill_ghosts_h_b_F77(f_real *hx, f_real *hy, f_real *hz,
			     f_real *ex, f_real *ey, f_real *ez,
			     f_real *jxi, f_real *jyi, f_real *jzi);


void INIT_basic_F77(void);

void PML_coeff_init_param_F77(void);

// ----------------------------------------------------------------------
// Wrappers to be called from C that call into Fortran

void
PSC_set_globals(struct psc *psc)
{
  struct psc_coeff *p = &psc->coeff;
  PSC_set_globals_F77(&p->cori, &p->alpha, &p->eta);
}

void
PSC_set_patch(struct psc *psc, int p)
{
  struct psc_patch *patch = &psc->patch[p];
  int imn[3], imx[3];
  for (int d = 0; d < 3; d++) {
    imn[d] = 0;
    imx[d] = patch->ldims[d] - 1;
  }
  PSC_set_patch_F77(imn, imx, psc->ibn, &psc->dt, patch->dx);
}

void
PSC_set_domain(struct psc *psc)
{
  struct psc_domain *p = &psc->domain;
  int imax[3];

  for (int d = 0; d < 3; d++) {
    imax[d] = p->gdims[d] - 1;
  }
  int use_pml_ = 0;
  int ilo[3] = {};

  // assert(p->bnd_part_lo==p->bnd_part_hi); FIX ME
  PSC_set_domain_F77(p->length, p->gdims, ilo, imax, p->bnd_fld_lo, p->bnd_fld_hi,
		     p->bnd_part_lo, psc->domain.np, psc->ibn, &use_pml_);
}

void
PSC_set_coeff(struct psc *psc)
{
  struct psc_coeff *p = &psc->coeff;
  PSC_set_coeff_F77(&p->beta,
		    &p->wl, &p->ld, &p->vos, &p->vt, &p->wp);
}

void
OUT_params_set(struct psc *psc)
{
  struct psc_coeff *p = &psc->coeff;
  OUT_params_set_F77(&p->np, &p->nnp);
}

void
SETUP_field()
{
  SETUP_field_F77();
}

static void
PSC_set_timestep(struct psc *psc)
{
  PSC_set_timestep_F77(&psc->timestep);
}

void
PSC_set_params(struct psc *psc)
{
  struct psc_param *p = &psc->prm;
  PSC_set_params_F77(&p->qq, &p->mm, &p->tt, &p->cc, &p->eps0,
		     &p->nmax, &p->lw, &p->i0, &p->n0, &p->e0, &p->b0,
		     &p->j0, &p->rho0, &p->phi0, &p->a0);
}

void
PIC_push_part_yz(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
#if 0
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_yz_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
		       flds[JXI], flds[JYI], flds[JZI],
		       flds[EX], flds[EY], flds[EZ],
		       flds[HX], flds[HY], flds[HZ]);
#else
  int *ib = pf->ib;
  PIC_push_part_yz_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
		       &F3_C(pf, JXI, ib[0], ib[1], ib[2]),
		       &F3_C(pf, JYI, ib[0], ib[1], ib[2]),
		       &F3_C(pf, JZI, ib[0], ib[1], ib[2]),
		       &F3_C(pf, EX , ib[0], ib[1], ib[2]),
		       &F3_C(pf, EY , ib[0], ib[1], ib[2]),
		       &F3_C(pf, EZ , ib[0], ib[1], ib[2]),
		       &F3_C(pf, HX , ib[0], ib[1], ib[2]),
		       &F3_C(pf, HY , ib[0], ib[1], ib[2]),
		       &F3_C(pf, HZ , ib[0], ib[1], ib[2]));
#endif
}

void
PIC_push_part_xy(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_xy_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
		       flds[JXI], flds[JYI], flds[JZI],
		       flds[EX], flds[EY], flds[EZ],
		       flds[HX], flds[HY], flds[HZ]);
}

void
PIC_push_part_xz(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_xz_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
		       flds[JXI], flds[JYI], flds[JZI],
		       flds[EX], flds[EY], flds[EZ],
		       flds[HX], flds[HY], flds[HZ]);
}

void
PIC_push_part_xyz(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_xyz_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
			flds[JXI], flds[JYI], flds[JZI],
			flds[EX], flds[EY], flds[EZ],
			flds[HX], flds[HY], flds[HZ]);
}

void
PIC_push_part_z(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_z_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
		      flds[JXI], flds[JYI], flds[JZI],
		      flds[EX], flds[EY], flds[EZ],
		      flds[HX], flds[HY], flds[HZ]);
}

void
PIC_push_part_z_vay(struct psc *psc, int p, struct psc_particles *prts, struct psc_fields *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  fields_fortran_t::real_t **flds = pf->data;
  PIC_push_part_z_vay_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1], &psc->p2A, &psc->p2B,
			  flds[JXI], flds[JYI], flds[JZI],
			  flds[EX], flds[EY], flds[EZ],
			  flds[HX], flds[HY], flds[HZ]);
}

void
PIC_randomize(struct psc_particles *prts)
{
  PIC_randomize_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1]);
}

void
PIC_bin_coll(struct psc_particles *prts)
{
  PIC_bin_coll_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1]);
}

void
PIC_find_cell_indices(struct psc_particles *prts)
{
  PIC_find_cell_indices_F77(&prts->n_part, &psc_particles_fortran(prts)->particles[-1]);
}

void
SET_param_pml(struct psc *psc)
{
  SET_param_pml_F77(&psc->pml.thick, &psc->pml.cushion, &psc->pml.size, &psc->pml.order);
}

void
INIT_basic()
{
  INIT_basic_F77();
}

void
PIC_msa_e(struct psc_fields *pf)
{
  assert(0); // FIXME, local fields are now 0-based
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_msa_e_F77(flds[EX], flds[EY], flds[EZ],
		flds[HX], flds[HY], flds[HZ],
		flds[JXI], flds[JYI], flds[JZI]);
}

void
PIC_msa_h(struct psc_fields *pf)
{
  assert(0); // FIXME, local fields are now 0-based
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_msa_h_F77(flds[EX], flds[EY], flds[EZ],
		flds[HX], flds[HY], flds[HZ],
		flds[JXI], flds[JYI], flds[JZI]);
}

void
PIC_pml_msa(struct psc_fields *pf)
{
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_pml_msa_F77(flds[EX], flds[EY], flds[EZ],
		  flds[HX], flds[HY], flds[HZ],
		  flds[DX], flds[DY], flds[DZ],
		  flds[BX], flds[BY], flds[BZ],
		  flds[JXI], flds[JYI], flds[JZI],
		  flds[EPS], flds[MU]);
}

void
PIC_msb_h(struct psc_fields *pf)
{
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_msb_h_F77(flds[EX], flds[EY], flds[EZ],
		flds[HX], flds[HY], flds[HZ],
		flds[JXI], flds[JYI], flds[JZI]);
}

void
PIC_msb_e(struct psc_fields *pf)
{
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_msb_e_F77(flds[EX], flds[EY], flds[EZ],
		flds[HX], flds[HY], flds[HZ],
		flds[JXI], flds[JYI], flds[JZI]);
}

void
PIC_pml_msb(struct psc_fields *pf)
{
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(ppsc);
  PIC_pml_msb_F77(flds[EX], flds[EY], flds[EZ],
		  flds[HX], flds[HY], flds[HZ],
		  flds[DX], flds[DY], flds[DZ],
		  flds[BX], flds[BY], flds[BZ],
		  flds[JXI], flds[JYI], flds[JZI],
		  flds[EPS], flds[MU]);
}

void
PIC_fill_ghosts_h_b(struct psc *psc, int p, struct psc_fields *pf)
{
  assert(0); // FIXME, local fields are now 0-based
  fields_fortran_t::real_t **flds = pf->data;
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_fill_ghosts_h_b_F77(flds[HX], flds[HY], flds[HZ],
			  flds[EX], flds[EY], flds[EZ],
			  flds[JXI], flds[JYI], flds[JZI]);
}

// ----------------------------------------------------------------------
// Wrappers to be called from Fortran that continue to C

f_real
C_p_pulse_x1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_x2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_x1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_x2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}


f_real
C_p_pulse_y1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_y2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_y1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_y2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}


f_real
C_p_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

void
PML_coeff_init_param(struct psc *psc)
{
  PML_coeff_init_param_F77();
}

void
psc_setup_fortran(struct psc *psc)
{
  INIT_basic();
  PSC_set_globals(psc);
  PSC_set_params(psc);
  PSC_set_coeff(psc);
  PSC_set_domain(psc);
  SET_param_pml(psc);
  PSC_set_patch(psc, 0);
  OUT_params_set(psc);
  PML_coeff_init_param(psc);
  SETUP_field();
}

