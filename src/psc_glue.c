
#include "psc.h"

#include <assert.h>

#define PSC_set_globals_F77 F77_FUNC_(psc_set_globals, PSC_SET_GLOBALS)
#define PSC_set_patch_F77 F77_FUNC_(psc_set_patch, PSC_SET_PATCH)
#define PIC_set_variables_F77 F77_FUNC(pic_set_variables,PIC_SET_VARIABLES)
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
#define INIT_grid_map_F77 F77_FUNC(init_grid_map,INIT_GRID_MAP)
#define CALC_densities_F77 F77_FUNC(calc_densities,CALC_DENSITIES)
#define FIELDS_alloc_F77 F77_FUNC_(fields_alloc, FIELDS_ALLOC)
#define FIELDS_free_F77 F77_FUNC_(fields_free, FIELDS_FREE)
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
#define PIC_msa_F77 F77_FUNC_(pic_msa, PIC_MSA)
#define PIC_msb_F77 F77_FUNC_(pic_msb, PIC_MSB)
#define PIC_pml_msa_F77 F77_FUNC_(pic_pml_msa, PIC_PML_MSA)
#define PIC_pml_msb_F77 F77_FUNC_(pic_pml_msb, PIC_PML_MSB)
#define SERV_read_1_F77 F77_FUNC_(serv_read_1, SERV_READ_1)
#define SERV_read_2_F77 F77_FUNC_(serv_read_2, SERV_READ_2)
#define SERV_write_F77 F77_FUNC_(serv_write, SERV_WRITE)
#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)

#define p_pulse_z1__F77 F77_FUNC(p_pulse_z1_,P_PULSE_Z1_)
#define s_pulse_z1__F77 F77_FUNC(s_pulse_z1_,S_PULSE_Z1_)
#define p_pulse_z2__F77 F77_FUNC(p_pulse_z2_,P_PULSE_Z2_)
#define s_pulse_z2__F77 F77_FUNC(s_pulse_z2_,S_PULSE_Z2_)

#define C_fields_alloc_cb_F77 F77_FUNC(c_fields_alloc_cb,C_FIELDS_ALLOC_CB)
#define C_p_pulse_z1_F77 F77_FUNC(c_p_pulse_z1,C_P_PULSE_Z1)
#define C_s_pulse_z1_F77 F77_FUNC(c_s_pulse_z1,C_S_PULSE_Z1)
#define C_p_pulse_z2_F77 F77_FUNC(c_p_pulse_z2,C_P_PULSE_Z2)
#define C_s_pulse_z2_F77 F77_FUNC(c_s_pulse_z2,C_S_PULSE_Z2)

void PSC_set_globals_F77(f_real *cori, f_real *alpha, f_real *eta);

void PSC_set_patch_F77(f_int *imn, f_int *imx, f_int *rd,
		       f_real *dt, f_real *dxyz);

void PIC_set_variables_F77(f_real *wl, f_real *wp, f_int *n);

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
void OUT_field_F77(void);
void OUT_part_F77(f_int *niloc, particle_fortran_t *p_niloc);
void SET_param_pml_F77(f_int *thick, f_int *cushion, f_int *size, f_int *order);
void INIT_grid_map_F77(void);
void CALC_densities_F77(f_int *niloc, particle_fortran_t *p_niloc,
			f_real *ne, f_real *ni, f_real *nn);
void FIELDS_alloc_F77(void);
void FIELDS_free_F77(void);

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
void PIC_msa_F77(f_real *ex, f_real *ey, f_real *ez,
		 f_real *hx, f_real *hy, f_real *hz,
		 f_real *jxi, f_real *jyi, f_real *jzi);
void PIC_msb_F77(f_real *ex, f_real *ey, f_real *ez,
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


void INIT_basic_F77(void);

f_real p_pulse_z1__F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt);
f_real s_pulse_z1__F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt);
f_real p_pulse_z2__F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt);
f_real s_pulse_z2__F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt);

// ----------------------------------------------------------------------
// Wrappers to be called from C that call into Fortran

void
PSC_set_globals()
{
  struct psc_coeff *p = &psc.coeff;
  PSC_set_globals_F77(&p->cori, &p->alpha, &p->eta);
}

void
PSC_set_patch(int p)
{
  struct psc_patch *patch = &psc.patch[p];
  int imx[3];
  for (int d = 0; d < 3; d++) {
    imx[d] = patch->off[d] + patch->ldims[d] - 1;
  }
  PSC_set_patch_F77(patch->off, imx, psc.ibn, &psc.dt, psc.dx);
}

static void
PIC_set_variables()
{
  PIC_set_variables_F77(&psc.coeff.wl, &psc.coeff.wp, &psc.timestep);
}

void
PIC_push_part_yz(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_yz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xy(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_xy_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xz(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_xz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xyz(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_xyz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_a(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_yz_a_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_b(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_yz_b_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_z(int patch, particles_fortran_t *pp, fields_fortran_t *pf)
{
  PSC_set_patch(patch);
  PIC_push_part_z_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
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
  PIC_set_variables();
  PIC_find_cell_indices_F77(&pp->n_part, &pp->particles[-1]);
}

void
OUT_field()
{
  OUT_field_F77();
}

void
OUT_part(particles_fortran_t *pp)
{
  OUT_part_F77(&pp->n_part, &pp->particles[-1]);
}

void
SET_param_pml()
{
  SET_param_pml_F77(&psc.pml.thick, &psc.pml.cushion, &psc.pml.size, &psc.pml.order);
}

real
PSC_p_pulse_z1(real xx, real yy, real zz, real tt)
{
  f_real _xx = xx, _yy = yy, _zz = zz, _tt = tt;
  return p_pulse_z1__F77(&_xx, &_yy, &_zz, &_tt);
}

real
PSC_s_pulse_z1(real xx, real yy, real zz, real tt)
{
  f_real _xx = xx, _yy = yy, _zz = zz, _tt = tt;
  return s_pulse_z1__F77(&_xx, &_yy, &_zz, &_tt);
}

real
PSC_p_pulse_z2(real xx, real yy, real zz, real tt)
{
  f_real _xx = xx, _yy = yy, _zz = zz, _tt = tt;
  return p_pulse_z2__F77(&_xx, &_yy, &_zz, &_tt);
}

real
PSC_s_pulse_z2(real xx, real yy, real zz, real tt)
{
  f_real _xx = xx, _yy = yy, _zz = zz, _tt = tt;
  return s_pulse_z2__F77(&_xx, &_yy, &_zz, &_tt);
}

void
CALC_densities(particles_fortran_t *pp, fields_fortran_t *pf)
{
  INIT_grid_map();
  SET_param_coeff();
  CALC_densities_F77(&pp->n_part, &pp->particles[-1],
		     pf->flds[NE], pf->flds[NI], pf->flds[NN]);
}

void
INIT_basic()
{
  INIT_basic_F77();
}

void
INIT_grid_map()
{
  INIT_basic();
  SET_param_domain();
  INIT_grid_map_F77();
}

void
PIC_fax(fields_fortran_t *pf, int m)
{
  INIT_grid_map();
  PIC_fax_F77(pf->flds[m]);
}

void
PIC_fay(fields_fortran_t *pf, int m)
{
  INIT_grid_map();
  PIC_fay_F77(pf->flds[m]);
}

void
PIC_faz(fields_fortran_t *pf, int m)
{
  INIT_grid_map();
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
  INIT_grid_map();
  f_int niloc_n;
  PIC_pex_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pex_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_pey(particles_fortran_t *pp)
{
  INIT_grid_map();
  f_int niloc_n;
  PIC_pey_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pey_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_pez(particles_fortran_t *pp)
{
  INIT_grid_map();
  f_int niloc_n;
  PIC_pez_a_F77(&pp->n_part, &pp->particles[-1], &niloc_n);
  particles_fortran_realloc(pp, niloc_n);
  PIC_pez_b_F77(&pp->n_part, &pp->particles[-1]);
}

void
PIC_msa(fields_fortran_t *pf)
{
  INIT_grid_map();
  SET_param_coeff();
  PIC_msa_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
	      pf->flds[HX], pf->flds[HY], pf->flds[HZ],
	      pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_pml_msa(fields_fortran_t *pf)
{
  INIT_grid_map();
  SET_param_coeff();
  PIC_pml_msa_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		  pf->flds[HX], pf->flds[HY], pf->flds[HZ],
		  pf->flds[DX], pf->flds[DY], pf->flds[DZ],
		  pf->flds[BX], pf->flds[BY], pf->flds[BZ],
		  pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		  pf->flds[EPS], pf->flds[MU]);
}

void
PIC_msb(fields_fortran_t *pf)
{
  PIC_msb_F77(pf->flds[EX], pf->flds[EY], pf->flds[EZ],
	      pf->flds[HX], pf->flds[HY], pf->flds[HZ],
	      pf->flds[JXI], pf->flds[JYI], pf->flds[JZI]);
}

void
PIC_pml_msb(fields_fortran_t *pf)
{
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
SERV_write(particles_fortran_t *pp, fields_fortran_t *pf)
{
  SERV_write_F77(&psc.timestep,
		 &pp->n_part, &pp->particles[-1],
		 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

// ----------------------------------------------------------------------
// Wrappers to be called from Fortran that continue to C

f_real
C_p_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  return psc_p_pulse_z1(*xx, *yy, *zz, *tt);
}


f_real
C_s_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  return psc_s_pulse_z1(*xx, *yy, *zz, *tt);
}

f_real
C_p_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  return psc_p_pulse_z2(*xx, *yy, *zz, *tt);
}


f_real
C_s_pulse_z2_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  return psc_s_pulse_z2(*xx, *yy, *zz, *tt);
}


// ----------------------------------------------------------------------
// same thing for allocating fields

f_real *__f_flds[NR_FIELDS];

f_real **
ALLOC_field()
{
  assert(psc.nr_patches == 1);
  SET_param_domain();
  PSC_set_patch(0);
  FIELDS_alloc_F77();
  // the callback function below will have magically been called,
  // setting __f_flds
  return __f_flds;
}

void
FREE_field(void)
{
  FIELDS_free_F77();
}

void
C_fields_alloc_cb_F77(f_real *ne, f_real *ni, f_real *nn,
		      f_real *jxi, f_real *jyi, f_real *jzi,
		      f_real *ex, f_real *ey, f_real *ez,
		      f_real *hx, f_real *hy, f_real *hz,
		      f_real *dx, f_real *dy, f_real *dz,
		      f_real *bx, f_real *by, f_real *bz,
		      f_real *eps, f_real *mu)
{
  __f_flds[NE] = ne;
  __f_flds[NI] = ni;
  __f_flds[NN] = nn;
  __f_flds[JXI] = jxi;
  __f_flds[JYI] = jyi;
  __f_flds[JZI] = jzi;
  __f_flds[EX] = ex;
  __f_flds[EY] = ey;
  __f_flds[EZ] = ez;
  __f_flds[HX] = hx;
  __f_flds[HY] = hy;
  __f_flds[HZ] = hz;
  __f_flds[DX] = dx;
  __f_flds[DY] = dy;
  __f_flds[DZ] = dz;
  __f_flds[BX] = bx;
  __f_flds[BY] = by;
  __f_flds[BZ] = bz;
  __f_flds[EPS] = eps;
  __f_flds[MU] = mu;
}
