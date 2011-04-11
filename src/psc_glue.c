
#include "psc.h"
#include "psc_glue.h"

#include "psc_push_fields.h"
#include "psc_bnd_fields_private.h"
#include "psc_pulse.h"

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
#define PIC_push_part_yz_a_F77 F77_FUNC(pic_push_part_yz_a,PIC_PUSH_PART_YZ_A)
#define PIC_push_part_yz_b_F77 F77_FUNC(pic_push_part_yz_b,PIC_PUSH_PART_YZ_B)
#define PIC_push_part_z_F77 F77_FUNC(pic_push_part_z,PIC_PUSH_PART_Z)
#define PIC_sort_F77 F77_FUNC(pic_sort,PIC_SORT)
#define PIC_randomize_F77 F77_FUNC(pic_randomize,PIC_RANDOMIZE)
#define PIC_bin_coll_F77 F77_FUNC(pic_bin_coll,PIC_BIN_COLL)
#define PIC_find_cell_indices_F77 F77_FUNC(pic_find_cell_indices,PIC_FIND_CELL_INDICES)
#define OUT_field_F77 F77_FUNC(out_field,OUT_FIELD)
#define OUT_part_F77 F77_FUNC(out_part,OUT_PART)
#define SET_param_pml_F77 F77_FUNC(set_param_pml,SET_PARAM_PML)
#define CALC_densities_F77 F77_FUNC(calc_densities,CALC_DENSITIES)
#define PIC_fax_F77 F77_FUNC(pic_fax, PIC_FAX)
#define PIC_fay_F77 F77_FUNC(pic_fay, PIC_FAY)
#define PIC_faz_F77 F77_FUNC(pic_faz, PIC_FAZ)
#define PIC_fex_F77 F77_FUNC(pic_fex, PIC_FEX)
#define PIC_fey_F77 F77_FUNC(pic_fey, PIC_FEY)
#define PIC_fez_F77 F77_FUNC(pic_fez, PIC_FEZ)
#define PIC_pex_a_F77 F77_FUNC(pic_pex_a, PIC_PEX_A)
#define PIC_pex_b_F77 F77_FUNC(pic_pex_b, PIC_PEX_B)
#define PIC_pey_a_F77 F77_FUNC(pic_pey_a, PIC_PEY_A)
#define PIC_pey_b_F77 F77_FUNC(pic_pey_b, PIC_PEY_B)
#define PIC_pez_a_F77 F77_FUNC(pic_pez_a, PIC_PEZ_A)
#define PIC_pez_b_F77 F77_FUNC(pic_pez_b, PIC_PEZ_B)
#define PIC_msa_e_F77 F77_FUNC_(pic_msa_e, PIC_MSA_E)
#define PIC_msa_h_F77 F77_FUNC_(pic_msa_h, PIC_MSA_H)
#define PIC_msb_h_F77 F77_FUNC_(pic_msb_h, PIC_MSB_H)
#define PIC_msb_e_F77 F77_FUNC_(pic_msb_e, PIC_MSB_E)
#define PIC_pml_msa_F77 F77_FUNC_(pic_pml_msa, PIC_PML_MSA)
#define PIC_pml_msb_F77 F77_FUNC_(pic_pml_msb, PIC_PML_MSB)
#define SERV_read_1_F77 F77_FUNC_(serv_read_1, SERV_READ_1)
#define SERV_read_2_F77 F77_FUNC_(serv_read_2, SERV_READ_2)
#define SERV_write_F77 F77_FUNC_(serv_write, SERV_WRITE)
#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)
#define PIC_fill_ghosts_h_b_F77 F77_FUNC_(pic_fill_ghosts_h_b, PIC_FILL_GHOSTS_H_B)

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

void PIC_push_part_yz_a_F77(f_int *niloc, particle_fortran_t *p_niloc,
			    f_real *p2A, f_real *p2B,
			    f_real *jxi, f_real *jyi, f_real *jzi,
			    f_real *ex, f_real *ey, f_real *ez,
			    f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_yz_b_F77(f_int *niloc, particle_fortran_t *p_niloc,
			    f_real *p2A, f_real *p2B,
			    f_real *jxi, f_real *jyi, f_real *jzi,
			    f_real *ex, f_real *ey, f_real *ez,
			    f_real *hx, f_real *hy, f_real *hz);

void PIC_push_part_z_F77(f_int *niloc, particle_fortran_t *p_niloc,
			 f_real *p2A, f_real *p2B,
			 f_real *jxi, f_real *jyi, f_real *jzi,
			 f_real *ex, f_real *ey, f_real *ez,
			 f_real *hx, f_real *hy, f_real *hz);

void PIC_sort_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_randomize_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_bin_coll_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_find_cell_indices_F77(f_int *niloc, particle_fortran_t *p_niloc);
void OUT_field_F77(f_real *ne, f_real *ni, f_real *nn,
		   f_real *jxi, f_real *jyi, f_real *jzi,
		   f_real *ex, f_real *ey, f_real *ez,
		   f_real *hx, f_real *hy, f_real *hz,
		   f_real *dvx, f_real *dvy, f_real *dvz,
		   f_real *bx, f_real *by, f_real *bz);
void OUT_part_F77(f_int *niloc, particle_fortran_t *p_niloc);
void SET_param_pml_F77(f_int *thick, f_int *cushion, f_int *size, f_int *order);
void CALC_densities_F77(f_int *niloc, particle_fortran_t *p_niloc,
			f_real *ne, f_real *ni, f_real *nn);

void PIC_fax_F77(f_real *f);
void PIC_fay_F77(f_real *f);
void PIC_faz_F77(f_real *f);
void PIC_fex_F77(f_real *f);
void PIC_fey_F77(f_real *f);
void PIC_fez_F77(f_real *f);
void PIC_pex_a_F77(f_int *niloc, particle_fortran_t *p_niloc, f_int *niloc_n);
void PIC_pex_b_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_pey_a_F77(f_int *niloc, particle_fortran_t *p_niloc, f_int *niloc_n);
void PIC_pey_b_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_pez_a_F77(f_int *niloc, particle_fortran_t *p_niloc, f_int *niloc_n);
void PIC_pez_b_F77(f_int *niloc, particle_fortran_t *p_niloc);
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
void SERV_read_1_F77(f_int *timestep, f_int *n_part);
void SERV_read_2_F77(f_int *niloc, particle_fortran_t *p_niloc,
		     f_real *jxi, f_real *jyi, f_real *jzi,
		     f_real *ex, f_real *ey, f_real *ez,
		     f_real *hx, f_real *hy, f_real *hz);
void SERV_write_F77(f_int *timestep, f_int *niloc, particle_fortran_t *p_niloc,
		    f_real *jxi, f_real *jyi, f_real *jzi,
		    f_real *ex, f_real *ey, f_real *ez,
		    f_real *hx, f_real *hy, f_real *hz);
void PIC_fill_ghosts_h_b_F77(f_real *hx, f_real *hy, f_real *hz,
			     f_real *ex, f_real *ey, f_real *ez,
			     f_real *jxi, f_real *jyi, f_real *jzi);


void INIT_basic_F77(void);

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
  int imx[3];
  for (int d = 0; d < 3; d++) {
    imx[d] = patch->off[d] + patch->ldims[d] - 1;
  }
  PSC_set_patch_F77(patch->off, imx, psc->ibn, &psc->dt, psc->dx);
}

void
PSC_set_domain(struct psc *psc)
{
  struct psc_domain *p = &psc->domain;
  int imax[3], np[3];

  mrc_domain_get_param_int3(psc->mrc_domain, "np", np);
  for (int d = 0; d < 3; d++) {
    imax[d] = p->gdims[d] - 1;
  }
  int use_pml_ = p->use_pml;
  int ilo[3] = {};
  PSC_set_domain_F77(p->length, p->gdims, ilo, imax, p->bnd_fld_lo, p->bnd_fld_hi,
		       p->bnd_part, np, psc->ibn, &use_pml_);
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
PIC_push_part_yz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_yz_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xy(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_xy_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_xz_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xyz(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_xyz_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
			pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_a(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_yz_a_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_b(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_yz_b_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_z(struct psc *psc, int p, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_push_part_z_F77(&pp->n_part, &pp->particles[-1], &psc->p2A, &psc->p2B,
		      pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		      pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		      pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_sort(particles_fortran_t *pp)
{
  PIC_sort_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_randomize(particles_fortran_t *pp)
{
  PIC_randomize_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_bin_coll(particles_fortran_t *pp)
{
  PIC_bin_coll_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_find_cell_indices(particles_fortran_t *pp)
{
  PIC_find_cell_indices_F77(&pp->n_part, &pp->particles[-1]);
}

void
OUT_field(fields_fortran_t *pf)
{
  OUT_field_F77(pf->flds[NE], pf->flds[NI], pf->flds[NN],
		pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		pf->flds[DX], pf->flds[DY], pf->flds[DZ],
		pf->flds[BX], pf->flds[BY], pf->flds[BZ]);
}

void
OUT_part(particles_fortran_t *pp)
{
  OUT_part_F77(&pp->n_part, &pp->particles[-1]);
}

void
SET_param_pml(struct psc *psc)
{
  SET_param_pml_F77(&psc->pml.thick, &psc->pml.cushion, &psc->pml.size, &psc->pml.order);
}

void
CALC_densities(particles_fortran_t *pp, fields_fortran_t *pf)
{
  CALC_densities_F77(&pp->n_part, &pp->particles[-1],
		     pf->flds[NE], pf->flds[NI], pf->flds[NN]);
}

void
INIT_basic()
{
  INIT_basic_F77();
}

void
PIC_fax(fields_fortran_t *pf, int m)
{
  PIC_fax_F77(pf->flds[m]);
}

void
PIC_fay(fields_fortran_t *pf, int m)
{
  PIC_fay_F77(pf->flds[m]);
}

void
PIC_faz(fields_fortran_t *pf, int m)
{
  PIC_faz_F77(pf->flds[m]);
}

void
PIC_fex(fields_fortran_t *pf, int m)
{
  PIC_fex_F77(pf->flds[m]);
}

void
PIC_fey(fields_fortran_t *pf, int m)
{
  PIC_fey_F77(pf->flds[m]);
}

void
PIC_fez(fields_fortran_t *pf, int m)
{
  PIC_fez_F77(pf->flds[m]);
}

void
PIC_pex(particles_fortran_t *pp)
{
  f_int niloc_n;
  PIC_pex_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pex_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_pey(particles_fortran_t *pp)
{
  f_int niloc_n;
  PIC_pey_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pey_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_pez(particles_fortran_t *pp)
{
  f_int niloc_n;
  PIC_pez_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pez_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_msa_e(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_msa_e_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_msa_h(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_msa_h_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_pml_msa(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_pml_msa_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		  pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		  pf->flds[DX], pf->flds[DY], pf->flds[DZ],
		  pf->flds[BX], pf->flds[BY], pf->flds[BZ],
		  pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		  pf->flds[EPS], pf->flds[MU]);
}

void
PIC_msb_h(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_msb_h_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_msb_e(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_msb_e_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_pml_msb(fields_fortran_t *pf)
{
  PSC_set_timestep(&psc);
  PIC_pml_msb_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		  pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		  pf->flds[DX], pf->flds[DY], pf->flds[DZ],
		  pf->flds[BX], pf->flds[BY], pf->flds[BZ],
		  pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		  pf->flds[EPS], pf->flds[MU]);
}

void
SERV_read_1(int *timestep, int *n_part)
{
  SERV_read_1_F77(timestep, n_part);
}

void
SERV_read_2(particles_fortran_t *pp, fields_fortran_t *pf)
{
  SERV_read_2_F77(&pp->n_part, &pp->particles[-1],
		  pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		  pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		  pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
SERV_write(struct psc *psc, particles_fortran_t *pp, fields_fortran_t *pf)
{
  SERV_write_F77(&psc->timestep,
		 &pp->n_part, &pp->particles[-1],
		 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_fill_ghosts_h_b(struct psc *psc, int p, fields_fortran_t *pf)
{
  PSC_set_timestep(psc);
  PSC_set_patch(psc, p);
  PIC_fill_ghosts_h_b_F77(pf->flds[HX], pf->flds[HY], pf->flds[HZ],
			  pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			  pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

// ----------------------------------------------------------------------
// Wrappers to be called from Fortran that continue to C

f_real
C_p_pulse_x1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_x2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_x1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_x2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_x2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}


f_real
C_p_pulse_y1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_y2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_y1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_y2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_y2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}


f_real
C_p_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z1;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z2;
  return psc_pulse_field_p(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z1;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
}

f_real
C_s_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(psc.push_fields);
  struct psc_pulse *pulse = bnd_fields->pulse_z2;
  return psc_pulse_field_s(pulse, *xx, *yy, *zz, *tt);
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
  SETUP_field();
}

