
#include "psc.h"

#include <assert.h>

#define PIC_set_variables_F77 F77_FUNC(pic_set_variables,PIC_SET_VARIABLES)
#define PIC_push_part_xz_F77 F77_FUNC(pic_push_part_xz,PIC_PUSH_PART_XZ)
#define PIC_push_part_yz_F77 F77_FUNC(pic_push_part_yz,PIC_PUSH_PART_YZ)
#define PIC_push_part_yz_a_F77 F77_FUNC(pic_push_part_yz_a,PIC_PUSH_PART_YZ_A)
#define PIC_push_part_yz_b_F77 F77_FUNC(pic_push_part_yz_b,PIC_PUSH_PART_YZ_B)
#define PIC_push_part_z_F77 F77_FUNC(pic_push_part_z,PIC_PUSH_PART_Z)
#define PIC_sort_1_F77 F77_FUNC(pic_sort_1,PIC_SORT_1)
#define PIC_randomize_F77 F77_FUNC(pic_randomize,PIC_RANDOMIZE)
#define PIC_bin_coll_F77 F77_FUNC(pic_bin_coll,PIC_BIN_COLL)
#define PIC_find_cell_indices_F77 F77_FUNC(pic_find_cell_indices,PIC_FIND_CELL_INDICES)
#define INIT_partition_F77 F77_FUNC_(init_partition, INIT_PARTITION)
#define INIT_idistr_F77 F77_FUNC_(init_idistr, INIT_IDISTR)
#define OUT_field_1_F77 F77_FUNC(out_field_1,OUT_FIELD_1)
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
#define ALLOC_field_F77 F77_FUNC_(alloc_field, ALLOC_FIELD)
#define FREE_particles_F77 F77_FUNC_(free_particles, FREE_PARTICLES)
#define FREE_field_F77 F77_FUNC_(free_field, FREE_FIELD)
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
#define INIT_basic_F77 F77_FUNC_(init_basic, INIT_BASIC)

#define p_pulse_z1__F77 F77_FUNC(p_pulse_z1_,P_PULSE_Z1_)

#define C_push_part_yz_F77 F77_FUNC(c_push_part_yz,C_PUSH_PART_YZ)
#define C_push_part_z_F77 F77_FUNC(c_push_part_z,C_PUSH_PART_Z)
#define C_sort_F77 F77_FUNC(c_sort,C_SORT)
#define C_out_field_F77 F77_FUNC(c_out_field,C_OUT_FIELD)
#define C_alloc_particles_cb_F77 F77_FUNC(c_alloc_particles_cb,C_ALLOC_PARTICLES_CB)
#define C_alloc_field_cb_F77 F77_FUNC(c_alloc_field_cb,C_ALLOC_FIELD_CB)
#define C_p_pulse_z1_F77 F77_FUNC(c_p_pulse_z1,C_P_PULSE_Z1)

void PIC_set_variables_F77(f_int *i1mn, f_int *i2mn, f_int *i3mn,
			   f_int *i1mx, f_int *i2mx, f_int *i3mx,
			   f_int *rd1, f_int *rd2, f_int *rd3,
			   f_int *i1n, f_int *i2n, f_int *i3n,
			   f_int *i1x, f_int *i2x, f_int *i3x,
			   f_real *cori, f_real *alpha, f_real *eta, 
			   f_real *dt, f_real *dx, f_real *dy, f_real *dz,
			   f_real *wl, f_real *wp, f_int *n);

void PIC_push_part_xz_F77(f_int *niloc, struct f_particle *p_niloc,
			  f_real *p2A, f_real *p2B,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *bx, f_real *by, f_real *bz);

void PIC_push_part_yz_F77(f_int *niloc, struct f_particle *p_niloc,
			  f_real *p2A, f_real *p2B,
			  f_real *ne, f_real *ni, f_real *nn,
			  f_real *jxi, f_real *jyi, f_real *jzi,
			  f_real *ex, f_real *ey, f_real *ez,
			  f_real *bx, f_real *by, f_real *bz);

void PIC_push_part_yz_a_F77(f_int *niloc, struct f_particle *p_niloc,
			    f_real *p2A, f_real *p2B,
			    f_real *ne, f_real *ni, f_real *nn,
			    f_real *jxi, f_real *jyi, f_real *jzi,
			    f_real *ex, f_real *ey, f_real *ez,
			    f_real *bx, f_real *by, f_real *bz);

void PIC_push_part_yz_b_F77(f_int *niloc, struct f_particle *p_niloc,
			    f_real *p2A, f_real *p2B,
			    f_real *ne, f_real *ni, f_real *nn,
			    f_real *jxi, f_real *jyi, f_real *jzi,
			    f_real *ex, f_real *ey, f_real *ez,
			    f_real *bx, f_real *by, f_real *bz);

void PIC_push_part_z_F77(f_int *niloc, struct f_particle *p_niloc,
			 f_real *p2A, f_real *p2B,
			 f_real *ne, f_real *ni, f_real *nn,
			 f_real *jxi, f_real *jyi, f_real *jzi,
			 f_real *ex, f_real *ey, f_real *ez,
			 f_real *bx, f_real *by, f_real *bz);

void PIC_sort_1_F77(f_int *niloc, struct f_particle *p_niloc);
void PIC_randomize_F77(f_int *niloc, struct f_particle *p_niloc);
void PIC_bin_coll_F77(f_int *niloc, struct f_particle *p_niloc);
void PIC_find_cell_indices_F77(f_int *niloc, struct f_particle *p_niloc);
void INIT_partition_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
			f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x,
			f_int *niloc_new);
void INIT_idistr_F77(f_int *part_label_off, f_int *rd1n, f_int *rd1x,
		     f_int *rd2n, f_int *rd2x, f_int *rd3n, f_int *rd3x);
void OUT_field_1_F77(void);
void OUT_part_F77(void);
void SET_param_pml_F77(f_int *thick, f_int *cushion, f_int *size, f_int *order);
void INIT_grid_map_F77(void);
void CALC_densities_F77(void);
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
void ALLOC_field_F77(void);
void FREE_particles_F77(void);
void FREE_field_F77(void);

void PIC_fax_F77(f_real *f);
void PIC_fay_F77(f_real *f);
void PIC_faz_F77(f_real *f);
void PIC_fex_F77(f_real *f);
void PIC_fey_F77(f_real *f);
void PIC_fez_F77(f_real *f);
void PIC_pex_F77(void);
void PIC_pey_F77(void);
void PIC_pez_F77(void);
void PIC_msa_F77(void);
void PIC_msb_F77(void);
void PIC_pml_msa_F77(void);
void PIC_pml_msb_F77(void);
void INIT_basic_F77(void);

f_real p_pulse_z1__F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt);

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
PIC_push_part_yz()
{
  PIC_set_variables();
  PIC_push_part_yz_F77(&psc.n_part, &psc.f_part[-1], &psc.p2A, &psc.p2B,
		       psc.f_fields[NE], psc.f_fields[NI], psc.f_fields[NN],
		       psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
		       psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
		       psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

void
PIC_push_part_xz()
{
  PIC_set_variables();
  PIC_push_part_xz_F77(&psc.n_part, &psc.f_part[-1], &psc.p2A, &psc.p2B,
		       psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
		       psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
		       psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

void
PIC_push_part_yz_a()
{
  PIC_set_variables();
  PIC_push_part_yz_a_F77(&psc.n_part, &psc.f_part[-1], &psc.p2A, &psc.p2B,
			 psc.f_fields[NE], psc.f_fields[NI], psc.f_fields[NN],
			 psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
			 psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
			 psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

void
PIC_push_part_yz_b()
{
  PIC_set_variables();
  PIC_push_part_yz_b_F77(&psc.n_part, &psc.f_part[-1], &psc.p2A, &psc.p2B,
			 psc.f_fields[NE], psc.f_fields[NI], psc.f_fields[NN],
			 psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
			 psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
			 psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

void
PIC_push_part_z()
{
  PIC_set_variables();
  PIC_push_part_z_F77(&psc.n_part, &psc.f_part[-1], &psc.p2A, &psc.p2B,
		      psc.f_fields[NE], psc.f_fields[NI], psc.f_fields[NN],
		      psc.f_fields[JXI], psc.f_fields[JYI], psc.f_fields[JZI],
		      psc.f_fields[EX], psc.f_fields[EY], psc.f_fields[EZ],
		      psc.f_fields[BX], psc.f_fields[BY], psc.f_fields[BZ]);
}

void
PIC_sort_1()
{
  PIC_sort_1_F77(&psc.n_part, &psc.f_part[-1]);
}

void
PIC_randomize()
{
  PIC_randomize_F77(&psc.n_part, &psc.f_part[-1]);
}

void
PIC_bin_coll()
{
  PIC_bin_coll_F77(&psc.n_part, &psc.f_part[-1]);
}

void
PIC_find_cell_indices()
{
  PIC_set_variables();
  PIC_find_cell_indices_F77(&psc.n_part, &psc.f_part[-1]);
}

void
OUT_field_1()
{
  OUT_field_1_F77();
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
CALC_densities()
{
  CALC_densities_F77();
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
PIC_fax(int m)
{
  INIT_grid_map();
  PIC_fax_F77(psc.f_fields[m]);
}

void
PIC_fay(int m)
{
  PIC_fay_F77(psc.f_fields[m]);
}

void
PIC_faz(int m)
{
  PIC_faz_F77(psc.f_fields[m]);
}

void
PIC_fex(int m)
{
  PIC_fex_F77(psc.f_fields[m]);
}

void
PIC_fey(int m)
{
  PIC_fey_F77(psc.f_fields[m]);
}

void
PIC_fez(int m)
{
  PIC_fez_F77(psc.f_fields[m]);
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
PIC_msa()
{
  INIT_grid_map();
  SET_param_coeff();
#ifdef USE_PML
  PIC_pml_msa_F77();
#else
  PIC_msa_F77();
#endif
}

void
PIC_msb()
{
#ifdef USE_PML
  PIC_pml_msb_F77(); 
#else
  PIC_msb_F77(); 
#endif
}

// ----------------------------------------------------------------------
// Wrappers to be called from Fortran that continue to C

void
C_push_part_yz_F77(f_real *p2A, f_real *p2B,
		   f_int *niloc, struct f_particle *p_niloc,
		   f_int *dummy)
{
  // make sure we got passed the right number of arguments
  assert(*dummy == 99);

  psc.p2A = *p2A;
  psc.p2B = *p2B;

  psc.n_part = *niloc;
  psc.f_part = &p_niloc[1];

  psc_push_part_yz();

  *p2A = psc.p2A;
  *p2A = psc.p2B;
}

void
C_push_part_z_F77(f_real *p2A, f_real *p2B,
		  f_int *niloc, struct f_particle *p_niloc,
		  f_int *dummy)
{
  // make sure we got passed the right number of arguments
  assert(*dummy == 99);

  psc.p2A = *p2A;
  psc.p2B = *p2B;

  psc.n_part = *niloc;
  psc.f_part = &p_niloc[1];

  psc_push_part_z();

  *p2A = psc.p2A;
  *p2A = psc.p2B;
}

void
C_sort_F77(f_int *niloc, struct f_particle *p_niloc,
	   f_int *dummy)
{
  // make sure we got passed the right number of arguments
  assert(*dummy == 99);

  psc.n_part = *niloc;
  psc.f_part = &p_niloc[1];

  psc_sort();
}

void
C_out_field_F77(f_int *dummy)
{
  // make sure we got passed the right number of arguments
  assert(*dummy == 99);

  psc_out_field();
}

f_real
C_p_pulse_z1_F77(f_real *xx, f_real *yy, f_real *zz, f_real *tt)
{
  return psc_p_pulse_z1(*xx, *yy, *zz, *tt);
}

// ----------------------------------------------------------------------
// slightly hacky way to do the equivalent of "malloc" using
// Fortran's allocate()

static struct f_particle *__f_part;

struct f_particle *
ALLOC_particles(int n_part)
{
  ALLOC_particles_F77(&n_part);
  // the callback function below will have magically been called,
  // setting __f_part
  return __f_part;
}

struct f_particle *
REALLOC_particles(int n_part_n)
{
  REALLOC_particles_F77(&n_part_n);
  // the callback function below will have magically been called,
  // setting __f_part
  return __f_part;
}

void
C_alloc_particles_cb_F77(struct f_particle *p_niloc)
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
  ALLOC_field_F77();
  // the callback function below will have magically been called,
  // setting __f_flds
  return __f_flds;
}

void
C_alloc_field_cb_F77(f_real *ne, f_real *ni, f_real *nn,
		     f_real *jxi, f_real *jyi, f_real *jzi,
		     f_real *ex, f_real *ey, f_real *ez,
		     f_real *bx, f_real *by, f_real *bz)
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
  __f_flds[BX] = bx;
  __f_flds[BY] = by;
  __f_flds[BZ] = bz;
}

void
FREE_particles(void)
{
  FREE_particles_F77();
}

void
FREE_field(void)
{
  FREE_field_F77();
}

