
#ifndef PSC_H
#define PSC_H

#include <config.h>

#include <mrc_domain.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#include "psc_pulse.h"

// ----------------------------------------------------------------------

#define BOUNDS_CHECK

// ----------------------------------------------------------------------

#define FIELDS_FORTRAN 1
#define FIELDS_C       2
#define FIELDS_SSE2    3

#define PARTICLES_FORTRAN 1
#define PARTICLES_C       2
#define PARTICLES_SSE2    3

// FIELDS_BASE and PARTICLES_BASE macros are defined by configure
// #define FIELDS_BASE FIELDS_FORTRAN
// #define PARTICLES_BASE PARTICLES_FORTRAN

enum {
  NE , NI , NN ,
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
  DX , DY , DZ ,
  BX , BY , BZ ,
  EPS, MU ,
  NR_FIELDS,
};

const char *fldname[NR_FIELDS];

// C floating point type
// used to switch between single and double precision

typedef float real;

#define real(x) x ## f

// Fortran types

typedef double f_real;
typedef int f_int;

// always need the fortran types for fortran interface
#include "psc_particles_fortran.h"
#include "psc_fields_fortran.h"

// ----------------------------------------------------------------------
// base particles type

#if PARTICLES_BASE == PARTICLES_FORTRAN

typedef particles_fortran_t particles_base_t;
typedef mparticles_fortran_t mparticles_base_t;
typedef particle_fortran_t particle_base_t;
typedef particle_fortran_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL MPI_PARTICLES_FORTRAN_REAL

#define particles_base_alloc   particles_fortran_alloc
#define particles_base_realloc particles_fortran_realloc
#define particles_base_free    particles_fortran_free
#define particles_base_get_one particles_fortran_get_one

#elif PARTICLES_BASE == PARTICLES_C

#include "psc_particles_c.h"

typedef particles_c_t particles_base_t;
typedef mparticles_c_t mparticles_base_t;
typedef particle_c_t particle_base_t;
typedef particle_c_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_C_REAL

#define particles_base_alloc   particles_c_alloc
#define particles_base_realloc particles_c_realloc
#define particles_base_free    particles_c_free
#define particles_base_get_one particles_c_get_one

#elif PARTICLES_BASE == PARTICLES_SSE2

#include "psc_particles_sse2.h"

typedef particles_sse2_t particles_base_t;
typedef mparticles_sse2_t mparticles_base_t;
typedef particle_sse2_t particle_base_t;
typedef particle_sse2_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_SSE2_REAL

#define particles_base_alloc   particles_sse2_alloc
#define particles_base_realloc particles_sse2_realloc
#define particles_base_free    particles_sse2_free
#define particles_base_get_one particles_sse2_get_one

#else
#error unknown PARTICLES_BASE
#endif

// ----------------------------------------------------------------------
// base fields type

#if FIELDS_BASE == FIELDS_FORTRAN

typedef fields_fortran_t fields_base_t;
typedef fields_fortran_real_t fields_base_real_t;
typedef mfields_fortran_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_FORTRAN_REAL

#define fields_base_alloc            fields_fortran_alloc
#define fields_base_alloc_with_array fields_fortran_alloc_with_array
#define fields_base_free             fields_fortran_free
#define fields_base_zero             fields_fortran_zero
#define fields_base_zero_all         fields_fortran_zero_all
#define fields_base_set              fields_fortran_set
#define fields_base_copy             fields_fortran_copy
#define fields_base_axpy_all         fields_fortran_axpy_all
#define fields_base_scale_all        fields_fortran_scale_all
#define fields_base_size             fields_fortran_size

#define F3_BASE(pf, m, jx,jy,jz)     F3_FORTRAN(pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_C

#include "psc_fields_c.h"

typedef fields_c_t fields_base_t;
typedef fields_c_real_t fields_base_real_t;
typedef mfields_c_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_C_REAL

#define fields_base_alloc            fields_c_alloc
#define fields_base_alloc_with_array fields_c_alloc_with_array
#define fields_base_free             fields_c_free
#define fields_base_zero             fields_c_zero
#define fields_base_zero_all         fields_c_zero_all
#define fields_base_set              fields_c_set
#define fields_base_copy             fields_c_copy
#define fields_base_axpy_all         fields_c_axpy_all
#define fields_base_scale_all        fields_c_scale_all
#define fields_base_size             fields_c_size

#define F3_BASE(pf, m, jx,jy,jz)     F3_C(pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_SSE2

#include "psc_fields_sse2.h"

typedef fields_sse2_t fields_base_t;
typedef fields_sse2_real_t fields_base_real_t;
typedef mfields_sse2_t mfields_base_t;
#define MPI_FIELDS_BASE_REAL MPI_FIELDS_SSE2_REAL

#define fields_base_alloc fields_sse2_alloc
#define fields_base_free  fields_sse2_free
#define fields_base_zero  fields_sse2_zero
#define fields_base_set   fields_sse2_set
#define fields_base_copy  fields_sse2_copy

#define F3_BASE(pf, m, jx,jy,jz) F3_SSE2(pf, m, jx,jy,jz)

#else
#error unknown FIELDS_BASE
#endif

#include "psc_case.h"

// user settable parameters
struct psc_param {
  double qq;
  double mm;
  double tt;
  double cc;
  double eps0;
  int nmax;
  double cpum;
  double lw;
  double i0;
  double n0;
  double e0;
  double b0;
  double j0;
  double rho0;
  double phi0;
  double a0;
  int nicell;
  bool const_num_particles_per_cell;
  bool fortran_particle_weight_hack;
  bool adjust_dt_to_cycles;
  double wallclock_limit;
  bool from_checkpoint;
};

// coefficients needed for computations
// -- derived, not provided by user
struct psc_coeff {
  double cori;
  double alpha;
  double beta;
  double eta;

  // FIXME are these needed in general?
  double wl;
  double ld;
  double vos;
  double vt;
  double wp;
  int np; // # steps for time-averaging fields
  int nnp; // # steps per laser cycle
};

struct psc_pml {
  int thick; // # grid points for PML
  int cushion; // # grid points for buffer zone
  int size; // # grid points PML + buffer
  int order; // PML order
};

// need to match fortran values

enum {
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_UPML,
  BND_FLD_TIME,
};

enum {
  BND_PART_REFLECTING,
  BND_PART_PERIODIC,
};

struct psc_domain {
  double length[3];
  int gdims[3];
  int bnd_fld_lo[3], bnd_fld_hi[3], bnd_part[3];
  bool use_pml;
};

void mfields_base_alloc(mfields_base_t *flds, int nr_fields);

// ----------------------------------------------------------------------
// general info / parameters for the code

struct psc_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*push_part_xy)(mfields_base_t *flds_base, mparticles_base_t *particles_base);
  void (*push_part_xz)(mfields_base_t *flds_base, mparticles_base_t *particles_base);
  void (*push_part_yz)(mfields_base_t *flds_base, mparticles_base_t *particles_base);
  void (*push_part_xyz)(mfields_base_t *flds_base, mparticles_base_t *particles_base);
  void (*push_part_z)(mfields_base_t *flds_base, mparticles_base_t *particles_base);
  void (*push_part_yz_a)(mfields_base_t *flds_base, mparticles_base_t *particles_base); // only does the simple first half step
  void (*push_part_yz_b)(mfields_base_t *flds_base, mparticles_base_t *particles_base); // 1/2 x and 1/1 p step
};

struct psc_push_field_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*push_field_a)(mfields_base_t *); // 1st half step
  void (*push_field_b)(mfields_base_t *); // 2nd half step
};

// FIXME, the randomize / sort interaction needs more work
// In particular, it's better to randomize just per-cell after the sorting

struct psc_randomize_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*randomize)(mparticles_base_t *);
};

struct psc_sort_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*sort)(mparticles_base_t *particles);
};

struct psc_collision_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*collision)(mparticles_base_t *particles);
};

struct psc_output_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*out_field)(mfields_base_t *flds, mparticles_base_t *particles);
  void (*dump_field)(int m, const char *fname);
  void (*dump_particles)(const char *fname);
};

struct psc_bnd_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*add_ghosts)(mfields_base_t *flds, int mb, int me);
  void (*fill_ghosts)(mfields_base_t *flds, int mb, int me);
  void (*exchange_particles)(mparticles_base_t *particles);
};

struct psc_moment_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*calc_densities)(mfields_base_t *pf_base, mparticles_base_t *pp_base,
			 mfields_base_t *pf);
  void (*calc_v)(mfields_base_t *pf_base, mparticles_base_t *pp_base,
		 mfields_base_t *pf);
  void (*calc_vv)(mfields_base_t *pf_base, mparticles_base_t *pp_base,
		  mfields_base_t *pf);
};

struct psc_patch {
  int ldims[3];       // size of local domain (w/o ghost points)
  int off[3];         // local to global offset
  double xb[3];       // lower left corner of the domain in this patch
};

#define CRDX(p, jx) (psc.dx[0] * ((jx) + psc.patch[p].off[0]))
#define CRDY(p, jy) (psc.dx[1] * ((jy) + psc.patch[p].off[1]))
#define CRDZ(p, jz) (psc.dx[2] * ((jz) + psc.patch[p].off[2]))

struct psc {
  struct psc_ops *ops;
  struct psc_push_field_ops *push_field_ops;
  struct psc_randomize_ops *randomize_ops;
  struct psc_sort_ops *sort_ops;
  struct psc_collision_ops *collision_ops;
  struct psc_output_ops *output_ops;
  struct psc_bnd_ops *bnd_ops;
  void *bnd_data;
  struct psc_moment_ops *moment_ops;
  struct psc_pulse *pulse_p_z1;
  struct psc_pulse *pulse_s_z1;
  struct psc_pulse *pulse_p_z2;
  struct psc_pulse *pulse_s_z2;
  struct psc_case *Case;
  // user-configurable parameters
  struct psc_param prm;
  struct psc_coeff coeff;
  struct psc_domain domain;
  struct psc_pml pml;

  // other parameters / constants
  double p2A, p2B;
  int timestep;
  double dt;
  double dx[3];

  mparticles_base_t particles;
  mfields_base_t flds;
  struct mrc_domain *mrc_domain;

  int nr_patches;
  struct psc_patch patch[1];
  int ibn[3];         // number of ghost points

  // C data structures
  void *c_ctx;

  // did we allocate the fields / particles (otherwise, Fortran did)
  bool allocated;

  double time_start;
};

#define foreach_3d(p, ix, iy, iz, l, r) {				\
  int __ilo[3] = { -l, -l, -l };					\
  int __ihi[3] = { psc.patch[p].ldims[0] + r,				\
		   psc.patch[p].ldims[1] + r,				\
		   psc.patch[p].ldims[2] + r };				\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define foreach_3d_end				\
  } } }

#define foreach_3d_g(p, ix, iy, iz) {					\
  int __ilo[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };		\
  int __ihi[3] = { psc.patch[p].ldims[0] + psc.ibn[0],			\
		   psc.patch[p].ldims[1] + psc.ibn[1],			\
		   psc.patch[p].ldims[2] + psc.ibn[2] };		\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define foreach_3d_g_end				\
  } } }

#define foreach_patch(p)				\
  for (int p = 0; p < psc.nr_patches; p++)


// ----------------------------------------------------------------------
// psc_config

struct psc_mod_config {
  const char *mod_particle;
  const char *mod_field;
  const char *mod_randomize;
  const char *mod_sort;
  const char *mod_collision;
  const char *mod_output;
  const char *mod_bnd;
  const char *mod_moment;
};

// we keep this info global for now.
// FIXME, I'd like to declare this extern, but mac os has a problem with that...

struct psc psc;

static inline void
psc_local_to_global_indices(int p, int jx, int jy, int jz,
			    int *ix, int *iy, int *iz)
{
  *ix = jx + psc.patch[p].off[0];
  *iy = jy + psc.patch[p].off[1];
  *iz = jz + psc.patch[p].off[2];
}

void psc_create(struct psc_mod_config *conf);
void psc_destroy(void);

void psc_init(const char *case_name);
void psc_init_param(const char *case_name);
void psc_init_partition(int *particle_label_offset);
void psc_init_particles(int particle_label_offset);
void psc_init_field(mfields_base_t *flds);
void psc_integrate(void);
void psc_push_particles(mfields_base_t *flds_base, mparticles_base_t *particles);
void psc_push_field_a(mfields_base_t *flds);
void psc_push_field_b(mfields_base_t *flds);
void psc_add_ghosts(mfields_base_t *flds, int mb, int me);
void psc_fill_ghosts(mfields_base_t *flds, int mb, int me);
void psc_exchange_particles(mparticles_base_t *particles);
void psc_calc_densities(mfields_base_t *flds, mparticles_base_t *particles,
			mfields_base_t *f);
void psc_calc_moments_v(mfields_base_t *flds, mparticles_base_t *particles,
			mfields_base_t *f);
void psc_calc_moments_vv(mfields_base_t *flds, mparticles_base_t *particles,
			 mfields_base_t *f);

void psc_dump_particles(mparticles_base_t *particles, const char *fname);
void psc_dump_field(mfields_base_t *flds, int m, const char *fname);

void psc_push_part_xyz(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_push_part_yz(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_push_part_z(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_push_part_xy(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_push_part_yz_a(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_push_part_yz_b(mfields_base_t *flds_base, mparticles_base_t *particles_base);
void psc_randomize(mparticles_base_t *particles);
void psc_sort(mparticles_base_t *particles);
void psc_collision(mparticles_base_t *particles);
void psc_out_field(mfields_base_t *flds, mparticles_base_t *particles);
void psc_out_particles(void);
void psc_set_n_particles(particles_base_t *pp, int n_part);

void psc_read_checkpoint(void);
void psc_write_checkpoint(void);

real psc_p_pulse_z1(real xx, real yy, real zz, real tt);
real psc_s_pulse_z1(real xx, real yy, real zz, real tt);

real psc_p_pulse_z2(real xx, real yy, real zz, real tt);
real psc_s_pulse_z2(real xx, real yy, real zz, real tt);


// various implementations of the psc
// (something like Fortran, generic C, CUDA, ...)

extern struct psc_ops psc_ops_fortran;
extern struct psc_ops psc_ops_generic_c;
extern struct psc_ops psc_ops_cuda;
extern struct psc_ops psc_ops_sse2; //Intel SIMD instructions
extern struct psc_ops psc_ops_none;

extern struct psc_push_field_ops psc_push_field_ops_fortran;
extern struct psc_push_field_ops psc_push_field_ops_c;
extern struct psc_push_field_ops psc_push_field_ops_none;

extern struct psc_randomize_ops psc_randomize_ops_fortran;
extern struct psc_randomize_ops psc_randomize_ops_none;

extern struct psc_sort_ops psc_sort_ops_fortran;
extern struct psc_sort_ops psc_sort_ops_qsort;
extern struct psc_sort_ops psc_sort_ops_countsort;
extern struct psc_sort_ops psc_sort_ops_countsort2;
extern struct psc_sort_ops psc_sort_ops_none;

extern struct psc_collision_ops psc_collision_ops_fortran;
extern struct psc_collision_ops psc_collision_ops_none;

extern struct psc_output_ops psc_output_ops_fortran;
extern struct psc_output_ops psc_output_ops_c;

extern struct psc_bnd_ops psc_bnd_ops_fortran;
extern struct psc_bnd_ops psc_bnd_ops_c;

extern struct psc_moment_ops psc_moment_ops_fortran;
extern struct psc_moment_ops psc_moment_ops_c;

extern struct psc_case_ops psc_case_ops_langmuir;
extern struct psc_case_ops psc_case_ops_wakefield;
extern struct psc_case_ops psc_case_ops_thinfoil;
extern struct psc_case_ops psc_case_ops_foils;
extern struct psc_case_ops psc_case_ops_curvedfoil;
extern struct psc_case_ops psc_case_ops_singlepart;
extern struct psc_case_ops psc_case_ops_harris;
extern struct psc_case_ops psc_case_ops_harris_xy;
extern struct psc_case_ops psc_case_ops_collisions;
extern struct psc_case_ops psc_case_ops_test_xz;
extern struct psc_case_ops psc_case_ops_test_yz;
extern struct psc_case_ops psc_case_ops_cone;

// Wrappers for Fortran functions
void PIC_push_part_xyz();
void PIC_push_part_xy(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_xz(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_z(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_a(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_push_part_yz_b(particles_fortran_t *pp, fields_fortran_t *pf);
void PIC_sort(particles_fortran_t *pp);
void PIC_randomize(particles_fortran_t *pp);
void PIC_bin_coll(particles_fortran_t *pp);
void PIC_find_cell_indices(particles_fortran_t *pp);
void PIC_msa(fields_fortran_t *pf);
void PIC_msb(fields_fortran_t *pf);
void PIC_pml_msa(fields_fortran_t *pf);
void PIC_pml_msb(fields_fortran_t *pf);
void OUT_field(void);
void OUT_part(void);
void CALC_densities(particles_fortran_t *pp, fields_fortran_t *pf);
void SET_param_domain(void);
void SET_param_pml(void);
void SET_param_psc(void);
void SET_param_coeff(void);
void SET_niloc(int niloc);
void SET_subdomain(void);
void GET_param_domain(void);
void GET_niloc(int *niloc);
void INIT_param_domain(void);
void INIT_param_psc(void);
void INIT_grid_map(void);
particle_fortran_t *ALLOC_particles(int n_part);
particle_fortran_t *REALLOC_particles(int n_part_n);
f_real **ALLOC_field(void);
void FREE_particles(void);
void FREE_field(void);
void INIT_basic(void);
void INIT_grid_map(void);
real PSC_p_pulse_z1(real x, real y, real z, real t);
real PSC_s_pulse_z1(real x, real y, real z, real t);
real PSC_p_pulse_z2(real x, real y, real z, real t);
real PSC_s_pulse_z2(real x, real y, real z, real t);

void PIC_fax(fields_fortran_t *pf, int m);
void PIC_fay(fields_fortran_t *pf, int m);
void PIC_faz(fields_fortran_t *pf, int m);
void PIC_fex(fields_fortran_t *pf, int m);
void PIC_fey(fields_fortran_t *pf, int m);
void PIC_fez(fields_fortran_t *pf, int m);
void PIC_pex(void);
void PIC_pey(void);
void PIC_pez(void);
void SERV_read_1(int *timestep, int *n_part);
void SERV_read_2(particles_fortran_t *pp, fields_fortran_t *pf);
void SERV_write(particles_fortran_t *pp, fields_fortran_t *pf);

// ----------------------------------------------------------------------
// other bits and hacks...

#define sqr(a) ((a) * (a))

#define HERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// ----------------------------------------------------------------------
// compiler bits

#ifndef __unused
#ifdef __GNUC__
#define	__unused	__attribute__((__unused__))
#else
#define	__unused	/* no attribute */
#endif
#endif

#endif
