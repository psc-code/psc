
#ifndef PSC_H
#define PSC_H

#include <config.h>

#include <stdbool.h>
#include <stdio.h>

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

// this matches the Fortran particle data structure

struct f_particle {
  f_real xi, yi, zi;
  f_real pxi, pyi, pzi;
  f_real qni;
  f_real mni;
  f_real cni;
  f_real lni;
  f_real wni;
};

// ----------------------------------------------------------------------
// macros to access Fortran fields

#define FF3_OFF(jx,jy,jz)						\
  (((((jz)-psc.ilg[2]))							\
    *psc.img[1] + ((jy)-psc.ilg[1]))					\
   *psc.img[0] + ((jx)-psc.ilg[0]))

#if 1

#define FF3(fldnr, jx,jy,jz)			\
  (psc.f_fields[fldnr][FF3_OFF(jx,jy,jz)])

#else

#define FF3(fldnr, jx,jy,jz)						\
  (*({int off = FF3_OFF(jx,jy,jz);					\
      assert(off >= 0);							\
      assert(off < psc.fld_size);					\
      &(psc.f_fields[fldnr][off]);					\
    }))

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
};

struct psc_domain {
  double length[3];
  int itot[3], ilo[3], ihi[3];
  int bnd_fld[3], bnd_part[3];
  int nghost[3];
  int nproc[3];
};

// ----------------------------------------------------------------------
// cases (initial conditions and other parameters...)

struct psc_case_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*init_param)(void);
  void (*init_field)(void);
  void (*init_nvt)(int kind, double x[3], double *q, double *m, double *n,
		   double v[3], double T[3]);
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
  void (*push_part_yz)(void);
  void (*push_part_z)(void);
  void (*push_part_yz_a)(void); // only does the simple first half step
  void (*push_part_yz_b)(void); // 1/2 x and 1/1 p step
};

struct psc_sort_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*sort)(void);
};

struct psc_output_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*out_field)(void);
};

struct psc {
  struct psc_ops *ops;
  struct psc_sort_ops *sort_ops;
  struct psc_output_ops *output_ops;
  struct psc_case_ops *case_ops;
  void *case_data;
  // user-configurable parameters
  struct psc_param prm;
  struct psc_coeff coeff;
  struct psc_domain domain;

  // other parameters / constants
  double p2A, p2B;
  int timestep;
  double dt;
  double dx[3];

  // Fortran compatible particles
  int n_part;
  struct f_particle *f_part;

  // Fortran compatible fields
  int ilo[3], ihi[3]; // local domain: il, il+1, ..., ih-1
  int ibn[3];         // number of ghost points
  int ilg[3], ihg[3]; // local domain incl ghost points: ilg, ilg+1, ..., ihg-1
  int img[3];         // total # points per dir incl. ghost points
  int fld_size;       // total # points per field incl. ghost points
  int glo[3], ghi[3]; // global domain

  f_real *f_fields[NR_FIELDS];

  // C data structures
  void *c_ctx;

  // did we allocate the fields / particles (otherwise, Fortran did)
  bool allocated;
};

// we keep this info global for now.
// FIXME, I'd like to declare this extern, but mac os has a problem with that...

struct psc psc;

void psc_create(const char *mod_particle, const char *mod_sort,
		const char *mod_output);
void psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part);
void psc_destroy();

void psc_init_param();
void psc_setup_parameters();
void psc_setup_fields_zero();
void psc_setup_fields_1();
void psc_setup_particles_1();
void psc_dump_particles(const char *fname);
void psc_save_particles_ref();
void psc_save_fields_ref();
void psc_check_currents_ref();
void psc_check_particles_ref();
void psc_check_particles_sorted();
void psc_create_test_1(const char *ops_name);

void psc_push_part_yz();
void psc_push_part_z();
void psc_push_part_yz_a();
void psc_push_part_yz_b();
void psc_sort();
void psc_out_field();

// various implementations of the psc
// (something like Fortran, generic C, CUDA, ...)

extern struct psc_ops psc_ops_fortran;
extern struct psc_ops psc_ops_generic_c;
extern struct psc_ops psc_ops_cuda;
extern struct psc_ops psc_ops_sse2; //Intel SIMD instructions

extern struct psc_sort_ops psc_sort_ops_fortran;
extern struct psc_sort_ops psc_sort_ops_qsort;
extern struct psc_sort_ops psc_sort_ops_countsort;
extern struct psc_sort_ops psc_sort_ops_countsort2;

extern struct psc_output_ops psc_output_ops_fortran;
extern struct psc_output_ops psc_output_ops_hdf5;

extern struct psc_case_ops psc_case_ops_harris;

// Wrappers for Fortran functions
void PIC_push_part_yz();
void PIC_push_part_z();
void PIC_push_part_yz_a();
void PIC_push_part_yz_b();
void PIC_sort_1();
void PIC_randomize();
void PIC_find_cell_indices();
void OUT_field_1();
void SET_param_domain();
void SET_param_psc();
void SET_param_coeff();
void INIT_param_domain();
void INIT_param_psc();

// ----------------------------------------------------------------------
// other bits and hacks...

#define sqr(a) ((a) * (a))

#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
