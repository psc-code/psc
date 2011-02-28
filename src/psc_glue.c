
#include "psc.h"

#include <assert.h>

#define PIC_set_variables_F77 F77_FUNC(pic_set_variables,PIC_SET_VARIABLES)
#define PIC_push_part_xy_F77 F77_FUNC(pic_push_part_xy,PIC_PUSH_PART_XY)
#define PIC_push_part_xz_F77 F77_FUNC(pic_push_part_xz,PIC_PUSH_PART_XZ)
#define PIC_push_part_yz_F77 F77_FUNC(pic_push_part_yz,PIC_PUSH_PART_YZ)
#define PIC_push_part_xyz_F77 F77_FUNC(pic_push_part_xyz,PIC_PUSH_PART_XYZ)
#define PIC_push_part_yz_a_F77 F77_FUNC(pic_push_part_yz_a,PIC_PUSH_PART_YZ_A)
#define PIC_push_part_yz_b_F77 F77_FUNC(pic_push_part_yz_b,PIC_PUSH_PART_YZ_B)
#define PIC_push_part_z_F77 F77_FUNC(pic_push_part_z,PIC_PUSH_PART_Z)
#define PIC_push_part_z_vay_F77 F77_FUNC(pic_push_part_z_vay,PIC_PUSH_PART_Z_VAY)
#define PIC_sort_F77 F77_FUNC(pic_sort,PIC_SORT)
#define PIC_randomize_F77 F77_FUNC(pic_randomize,PIC_RANDOMIZE)
#define PIC_bin_coll_F77 F77_FUNC(pic_bin_coll,PIC_BIN_COLL)
#define PIC_find_cell_indices_F77 F77_FUNC(pic_find_cell_indices,PIC_FIND_CELL_INDICES)
#define INIT_partition_F77 F77_FUNC_(init_partition, INIT_PARTITION)
#define INIT_idistr_F77 F77_FUNC_(init_idistr, INIT_IDISTR)
#define OUT_field_F77 F77_FUNC(out_field,OUT_FIELD)
#define OUT_part_F77 F77_FUNC(out_part,OUT_PART)
#define SET_param_pml_F77 F77_FUNC(set_param_pml,SET_PARAM_PML)
#define INIT_grid_map_F77 F77_FUNC(init_grid_map,INIT_GRID_MAP)
#define CALC_densities_F77 F77_FUNC(calc_densities,CALC_DENSITIES)
#define SET_subdomain_F77 F77_FUNC_(set_subdomain, SET_SUBDOMAIN)
#define GET_subdomain_F77 F77_FUNC_(get_subdomain, GET_SUBDOMAIN)
#define SET_niloc_F77 F77_FUNC_(set_niloc, SET_NILOC)
#define GET_niloc_F77 F77_FUNC_(get_niloc, GET_NILOC)
#define ALLOC_particles_F77 F77_FUNC_(alloc_particles, ALLOC_PARTICLES)
#define REALLOC_particles_F77 F77_FUNC_(realloc_particles, REALLOC_PARTICLES)
#define FREE_particles_F77 F77_FUNC_(free_particles, FREE_PARTICLES)
#define FIELDS_alloc_F77 F77_FUNC_(fields_alloc, FIELDS_ALLOC)
#define FIELDS_free_F77 F77_FUNC_(fields_free, FIELDS_FREE)
#define PIC_fax_F77 F77_FUNC(pic_fax, PIC_FAX)
#define PIC_fay_F77 F77_FUNC(pic_fay, PIC_FAY)
#define PIC_faz_F77 F77_FUNC(pic_faz, PIC_FAZ)
#define PIC_fex_F77 F77_FUNC(pic_fex, PIC_FEX)
#define PIC_fey_F77 F77_FUNC(pic_fey, PIC_FEY)
#define PIC_fez_F77 F77_FUNC(pic_fez, PIC_FEZ)
#define PIC_pex_F77 F77_FUNC(pic_pex, PIC_PEX)
#define PIC_pey_F77 F77_FUNC(pic_pey, PIC_PEY)
#define PIC_pez_F77 F77_FUNC(pic_pez, PIC_PEZ)
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

#define C_alloc_particles_cb_F77 F77_FUNC(c_alloc_particles_cb,C_ALLOC_PARTICLES_CB)
#define C_fields_alloc_cb_F77 F77_FUNC(c_fields_alloc_cb,C_FIELDS_ALLOC_CB)
#define C_p_pulse_z1_F77 F77_FUNC(c_p_pulse_z1,C_P_PULSE_Z1)
#define C_s_pulse_z1_F77 F77_FUNC(c_s_pulse_z1,C_S_PULSE_Z1)
#define C_p_pulse_z2_F77 F77_FUNC(c_p_pulse_z2,C_P_PULSE_Z2)
#define C_s_pulse_z2_F77 F77_FUNC(c_s_pulse_z2,C_S_PULSE_Z2)

void PIC_set_variables_F77(f_int *i1mn, f_int *i2mn, f_int *i3mn,
			   f_int *i1mx, f_int *i2mx, f_int *i3mx,
			   f_int *rd1, f_int *rd2, f_int *rd3,
			   f_int *i1n, f_int *i2n, f_int *i3n,
			   f_int *i1x, f_int *i2x, f_int *i3x,
			   f_real *cori, f_real *alpha, f_real *eta, 
			   f_real *dt, f_real *dx, f_real *dy, f_real *dz,
			   f_real *wl, f_real *wp, f_int *n);

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

void PIC_push_part_z_vay_F77(f_int *niloc, particle_fortran_t *p_niloc,
												 f_real *p2A, f_real *p2B,
												 f_real *jxi, f_real *jyi, f_real *jzi,
												 f_real *ex, f_real *ey, f_real *ez,
												 f_real *hx, f_real *hy, f_real *hz);

void PIC_sort_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_randomize_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_bin_coll_F77(f_int *niloc, particle_fortran_t *p_niloc);
void PIC_find_cell_indices_F77(f_int *niloc, particle_fortran_t *p_niloc);
void INIT_partition_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
			f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x,
			f_int *niloc_new);
void INIT_idistr_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
		     f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x);
void OUT_field_F77(void);
void OUT_part_F77(void);
void SET_param_pml_F77(f_int *thick, f_int *cushion, f_int *size, f_int *order);
void INIT_grid_map_F77(void);
void CALC_densities_F77(f_int *niloc, particle_fortran_t *p_niloc,
			f_real *ne, f_real *ni, f_real *nn);
void SET_subdomain_F77(f_int *i1mn, f_int *i1mx, f_int *i2mn, f_int *i2mx,
		       f_int *i3mn, f_int *i3mx, f_int *i1bn, f_int *i2bn,
		       f_int *i3bn);
void GET_subdomain_F77(f_int *i1mn, f_int *i1mx, f_int *i2mn, f_int *i2mx,
		       f_int *i3mn, f_int *i3mx, f_int *i1bn, f_int *i2bn,
		       f_int *i3bn);
void SET_niloc_F77(f_int *niloc);
void GET_niloc_F77(f_int *niloc);
void ALLOC_particles_F77(f_int *n_part);
void REALLOC_particles_F77(f_int *n_part_n);
void FREE_particles_F77(void);
void FIELDS_alloc_F77(void);
void FIELDS_free_F77(void);

void PIC_fax_F77(f_real *f);
void PIC_fay_F77(f_real *f);
void PIC_faz_F77(f_real *f);
void PIC_fex_F77(f_real *f);
void PIC_fey_F77(f_real *f);
void PIC_fez_F77(f_real *f);
void PIC_pex_F77(void);
void PIC_pey_F77(void);
void PIC_pez_F77(void);
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

static void
PIC_set_variables()
{
  int i0mx = psc.ihi[0] - 1, i1mx = psc.ihi[1] - 1, i2mx = psc.ihi[2] - 1;
  int i0x = psc.domain.ihi[0] - 1;
  int i1x = psc.domain.ihi[1] - 1;
  int i2x = psc.domain.ihi[2] - 1;

  PIC_set_variables_F77(&psc.ilo[0], &psc.ilo[1], &psc.ilo[2],
			&i0mx, &i1mx, &i2mx,
			&psc.ibn[0], &psc.ibn[1], &psc.ibn[2],
			&psc.domain.ilo[0], &psc.domain.ilo[1], &psc.domain.ilo[2],
			&i0x, &i1x, &i2x,
			&psc.coeff.cori, &psc.coeff.alpha, &psc.coeff.eta,
			&psc.dt, &psc.dx[0], &psc.dx[1], &psc.dx[2],
			&psc.coeff.wl, &psc.coeff.wp, &psc.timestep);
}

void
PIC_push_part_yz(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_yz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xy(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_xy_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xz(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_xz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		       pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		       pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		       pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_xyz(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_xyz_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_a(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_yz_a_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_yz_b(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_yz_b_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
			 pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
			 pf->flds[EX], pf->flds[EY], pf->flds[EZ],
			 pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}

void
PIC_push_part_z(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_z_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
		      pf->flds[JXI], pf->flds[JYI], pf->flds[JZI],
		      pf->flds[EX], pf->flds[EY], pf->flds[EZ],
		      pf->flds[HX], pf->flds[HY], pf->flds[HZ]);
}


void
PIC_push_part_z_vay(particles_fortran_t *pp, fields_fortran_t *pf)
{
  PIC_set_variables();
  PIC_push_part_z_vay_F77(&pp->n_part, &pp->particles[-1], &psc.p2A, &psc.p2B,
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
OUT_part()
{
  OUT_part_F77();
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


static int rd1n, rd1x, rd2n, rd2x, rd3n, rd3x;
static int part_label_offset;

void
INIT_partition(int *n_part)
{
  INIT_partition_F77(&part_label_offset, &rd1n, &rd1x, &rd2n, &rd2x, &rd3n, &rd3x,
		     n_part);
  GET_subdomain();
}

void
INIT_idistr(void)
{
  INIT_idistr_F77(&part_label_offset, &rd1n, &rd1x, &rd2n, &rd2x, &rd3n, &rd3x);
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
SET_subdomain()
{
  f_int i1mn = psc.ilo[0];
  f_int i2mn = psc.ilo[1];
  f_int i3mn = psc.ilo[2];
  f_int i1mx = psc.ihi[0] - 1;
  f_int i2mx = psc.ihi[1] - 1;
  f_int i3mx = psc.ihi[2] - 1;
  f_int i1bn = psc.ibn[0];
  f_int i2bn = psc.ibn[1];
  f_int i3bn = psc.ibn[2];
  SET_subdomain_F77(&i1mn, &i1mx, &i2mn, &i2mx, &i3mn, &i3mx, &i1bn, &i2bn, &i3bn);
}

void
GET_subdomain()
{
  f_int i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, i1bn, i2bn, i3bn;
  GET_subdomain_F77(&i1mn, &i1mx, &i2mn, &i2mx, &i3mn, &i3mx, &i1bn, &i2bn, &i3bn);
  psc.ilo[0] = i1mn;
  psc.ilo[1] = i2mn;
  psc.ilo[2] = i3mn;
  psc.ihi[0] = i1mx + 1;
  psc.ihi[1] = i2mx + 1;
  psc.ihi[2] = i3mx + 1;
  psc.ibn[0] = i1bn;
  psc.ibn[1] = i2bn;
  psc.ibn[2] = i3bn;
  for (int d = 0; d < 3; d++) {
    psc.ilg[d] = psc.ilo[d] - psc.ibn[d];
    psc.ihg[d] = psc.ihi[d] + psc.ibn[d];
    psc.img[d] = psc.ihg[d] - psc.ilg[d];
  }
  psc.fld_size = psc.img[0] * psc.img[1] * psc.img[2];
}

void
SET_niloc(int niloc)
{
  SET_niloc_F77(&niloc);
}

void
GET_niloc(int *p_niloc)
{
  GET_niloc_F77(p_niloc);
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
PIC_pex()
{
  PIC_pex_F77();
}

void
PIC_pey()
{
  PIC_pey_F77();
}

void
PIC_pez()
{
  PIC_pez_F77();
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
// slightly hacky way to do the equivalent of "malloc" using
// Fortran's allocate()

static particle_fortran_t *__f_part;

particle_fortran_t *
ALLOC_particles(int n_part)
{
  ALLOC_particles_F77(&n_part);
  // the callback function below will have magically been called,
  // setting __f_part
  return __f_part;
}

particle_fortran_t *
REALLOC_particles(int n_part_n)
{
  REALLOC_particles_F77(&n_part_n);
  // the callback function below will have magically been called,
  // setting __f_part
  return __f_part;
}

void
C_alloc_particles_cb_F77(particle_fortran_t *p_niloc)
{
  __f_part = &p_niloc[1];
}

// ----------------------------------------------------------------------
// same thing for allocating fields

f_real *__f_flds[NR_FIELDS];

f_real **
ALLOC_field()
{
  SET_param_domain();
  SET_subdomain();
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

void
FREE_particles(void)
{
  FREE_particles_F77();
}

