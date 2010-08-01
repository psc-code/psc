
#ifndef PSC_H
#define PSC_H

#include <config.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>

#include "psc_pulse.h"
#include "psc_case.h"

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
  BX , BY , BZ ,
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

// ----------------------------------------------------------------------
// macros to access Fortran fields

#define FF3_OFF(jx,jy,jz)						\
  (((((jz)-psc.ilg[2]))							\
    *psc.img[1] + ((jy)-psc.ilg[1]))					\
   *psc.img[0] + ((jx)-psc.ilg[0]))

#define _FF3(fld, jx,jy,jz)  (fld[FF3_OFF(jx,jy,jz)])

// always need the fortran types for fortran interface
#include "psc_particles_fortran.h"
#include "psc_fields_fortran.h"

// ----------------------------------------------------------------------
// base particles type

#if PARTICLES_BASE == PARTICLES_FORTRAN

typedef psc_particles_fortran_t psc_particles_base_t;
typedef particle_fortran_t particle_base_t;
typedef particle_fortran_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL MPI_PARTICLES_FORTRAN_REAL

#define psc_particles_base_alloc   psc_particles_fortran_alloc
#define psc_particles_base_realloc psc_particles_fortran_realloc
#define psc_particles_base_free    psc_particles_fortran_free
#define psc_particles_base_get_one psc_particles_fortran_get_one

#elif PARTICLES_BASE == PARTICLES_C

#include "psc_particles_c.h"

typedef psc_particles_c_t psc_particles_base_t;
typedef particle_c_t particle_base_t;
typedef particle_c_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_C_REAL

#define psc_particles_base_alloc   psc_particles_c_alloc
#define psc_particles_base_realloc psc_particles_c_realloc
#define psc_particles_base_free    psc_particles_c_free
#define psc_particles_base_get_one psc_particles_c_get_one

#elif PARTICLES_BASE == PARTICLES_SSE2

#include "psc_particles_sse2.h"

typedef psc_particles_sse2_t psc_particles_base_t;
typedef particle_sse2_t particle_base_t;
typedef particle_sse2_real_t particle_base_real_t;
#define MPI_PARTICLES_BASE_REAL    MPI_PARTICLES_SSE2_REAL

#define psc_particles_base_alloc   psc_particles_sse2_alloc
#define psc_particles_base_realloc psc_particles_sse2_realloc
#define psc_particles_base_free    psc_particles_sse2_free
#define psc_particles_base_get_one psc_particles_sse2_get_one

#else
#error unknown PARTICLES_BASE
#endif

// ----------------------------------------------------------------------
// base fields type

#if FIELDS_BASE == FIELDS_FORTRAN

typedef psc_fields_fortran_t psc_fields_base_t;
typedef fields_fortran_real_t fields_base_real_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_FORTRAN_REAL

#define psc_fields_base_alloc  psc_fields_fortran_alloc
#define psc_fields_base_free   psc_fields_fortran_free
#define psc_fields_base_zero   psc_fields_fortran_zero

#define F3_BASE(m, jx,jy,jz)  F3_FORTRAN(&psc.pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_C

#include "psc_fields_c.h"

typedef psc_fields_c_t psc_fields_base_t;
typedef fields_c_real_t fields_base_real_t;
#define MPI_FIELDS_BASE_REAL  MPI_FIELDS_C_REAL

#define psc_fields_base_alloc psc_fields_c_alloc
#define psc_fields_base_free  psc_fields_c_free
#define psc_fields_base_zero  psc_fields_c_zero

#define F3_BASE(m, jx,jy,jz)  F3_C(&psc.pf, m, jx,jy,jz)

#elif FIELDS_BASE == FIELDS_SSE2

#include "psc_fields_sse2.h"

typedef psc_fields_sse2_t psc_fields_base_t;
typedef fields_sse2_real_t fields_base_real_t;
#define MPI_FIELDS_BASE_REAL MPI_FIELDS_SSE2_REAL

#define psc_fields_base_alloc psc_fields_sse2_alloc
#define psc_fields_base_free  psc_fields_sse2_free
#define psc_fields_base_zero  psc_fields_sse2_zero

#define F3_BASE(m, jx,jy,jz)  F3_SSE2(&psc.pf, m, jx,jy,jz)

#else
#error unknown FIELDS_BASE
#endif

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
  int itot[3], ilo[3], ihi[3];
  int bnd_fld_lo[3], bnd_fld_hi[3], bnd_part[3];
  int nghost[3];
  int nproc[3];
};

// ----------------------------------------------------------------------
// general info / parameters for the code

struct psc_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*particles_from_fortran)(void);
  void (*particles_to_fortran)(void);
  void (*fields_from_fortran)(void);
  void (*fields_to_fortran)(void);
  void (*push_part_xz)(void);
  void (*push_part_yz)(void);
  void (*push_part_z)(void);
  void (*push_part_yz_a)(void); // only does the simple first half step
  void (*push_part_yz_b)(void); // 1/2 x and 1/1 p step
};

struct psc_push_field_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*push_field_a)(void); // 1st half step
  void (*push_field_b)(void); // 2nd half step
};

// FIXME, the randomize / sort interaction needs more work
// In particular, it's better to randomize just per-cell after the sorting

struct psc_randomize_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*randomize)(void);
};

struct psc_sort_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*sort)(void);
};

struct psc_collision_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*collision)(void);
};

struct psc_output_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*out_field)(void);
  void (*dump_field)(int m, const char *fname);
  void (*dump_particles)(const char *fname);
};

struct psc_bnd_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*add_ghosts)(int mb, int me);
  void (*fill_ghosts)(int mb, int me);
  void (*exchange_particles)(void);
};

struct psc_moment_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*calc_densities)(void);
};

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

  psc_particles_base_t pp;
  psc_fields_base_t pf;

  // Fortran compatible fields
  int ilo[3], ihi[3]; // local domain: il, il+1, ..., ih-1
  int ibn[3];         // number of ghost points
  int ilg[3], ihg[3]; // local domain incl ghost points: ilg, ilg+1, ..., ihg-1
  int img[3];         // total # points per dir incl. ghost points
  int fld_size;       // total # points per field incl. ghost points

  // C data structures
  void *c_ctx;

  // did we allocate the fields / particles (otherwise, Fortran did)
  bool allocated;
};

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

void psc_create(struct psc_mod_config *conf);
void psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part);
void psc_destroy();

void psc_init(const char *case_name);
void psc_init_param(const char *case_name);
void psc_init_partition(int *n_part, int *particle_label_offset);
void psc_init_particles(int particle_label_offset);
void psc_init_field();
void psc_integrate();
void psc_push_particles();
void psc_push_field_a();
void psc_push_field_b();
void psc_add_ghosts(int mb, int me);
void psc_fill_ghosts(int mb, int me);
void psc_exchange_particles(void);
void psc_calc_densities(void);

void psc_setup_parameters();
void psc_setup_fields_zero();
void psc_setup_fields_1();
void psc_setup_particles_1();
void psc_dump_particles(const char *fname);
void psc_dump_field(int m, const char *fname);
void psc_save_particles_ref();
void psc_save_fields_ref();
void psc_check_currents_ref(double thres);
void psc_check_currents_ref_noghost(double thres);
void psc_check_fields_ref(int *flds, double thres);
void psc_check_particles_ref(double thres);
void psc_check_particles_sorted();
void psc_create_test_1(const char *ops_name);
void psc_create_test_xz(struct psc_mod_config *conf);
void psc_create_test_yz(struct psc_mod_config *conf);

void psc_push_part_yz();
void psc_push_part_z();
void psc_push_part_yz_a();
void psc_push_part_yz_b();
void psc_randomize();
void psc_sort();
void psc_collision();
void psc_out_field();
void psc_out_particles();
void psc_set_n_particles(int n_part);

real psc_p_pulse_z1(real xx, real yy, real zz, real tt);

// various implementations of the psc
// (something like Fortran, generic C, CUDA, ...)

extern struct psc_ops psc_ops_fortran;
extern struct psc_ops psc_ops_generic_c;
extern struct psc_ops psc_ops_cuda;
extern struct psc_ops psc_ops_sse2; //Intel SIMD instructions

extern struct psc_push_field_ops psc_push_field_ops_fortran;
extern struct psc_push_field_ops psc_push_field_ops_c;

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
extern struct psc_moment_ops psc_moment_ops_generic_c;

extern struct psc_case_ops psc_case_ops_langmuir;
extern struct psc_case_ops psc_case_ops_wakefield;
extern struct psc_case_ops psc_case_ops_thinfoil;
extern struct psc_case_ops psc_case_ops_singlepart;
extern struct psc_case_ops psc_case_ops_harris;
extern struct psc_case_ops psc_case_ops_test_xz;
extern struct psc_case_ops psc_case_ops_test_yz;

// Wrappers for Fortran functions
void PIC_push_part_xz();
void PIC_push_part_yz();
void PIC_push_part_z();
void PIC_push_part_yz_a();
void PIC_push_part_yz_b();
void PIC_sort();
void PIC_randomize();
void PIC_bin_coll();
void PIC_find_cell_indices();
void PIC_msa();
void PIC_msb();
void OUT_field();
void OUT_part();
void CALC_densities(psc_fields_fortran_t *pf);
void SET_param_domain();
void SET_param_pml();
void SET_param_psc();
void SET_param_coeff();
void SET_niloc(int niloc);
void SET_subdomain();
void GET_param_domain();
void GET_niloc(int *niloc);
void GET_subdomain();
void INIT_param_domain();
void INIT_param_psc();
void INIT_grid_map();
void INIT_partition(int *n_part);
void INIT_idistr();
particle_fortran_t *ALLOC_particles(int n_part);
particle_fortran_t *REALLOC_particles(int n_part_n);
f_real **ALLOC_field();
void FREE_particles(void);
void FREE_field(void);
void INIT_basic(void);
void INIT_grid_map(void);
real PSC_p_pulse_z1(real x, real y, real z, real t);

void PIC_fax(psc_fields_fortran_t *pf, int m);
void PIC_fay(psc_fields_fortran_t *pf, int m);
void PIC_faz(psc_fields_fortran_t *pf, int m);
void PIC_fex(psc_fields_fortran_t *pf, int m);
void PIC_fey(psc_fields_fortran_t *pf, int m);
void PIC_fez(psc_fields_fortran_t *pf, int m);
void PIC_pex(void);
void PIC_pey(void);
void PIC_pez(void);
void PIC_msa(void);
void PIC_msb(void);

// ----------------------------------------------------------------------
// other bits and hacks...

#define sqr(a) ((a) * (a))

#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
