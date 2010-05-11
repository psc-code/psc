
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

// C floating point type
// used to switch between single and double precision

typedef float real;

// Fortran types

typedef double f_real;
typedef int f_int;

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

struct psc_param {
  double cori, eta, alpha;
  double wl;
  double wp;
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
  void (*push_part_yz_b)(void); // does everything but the currents
};

struct psc {
  struct psc_ops *ops;
  // user-configurable parameters
  struct psc_param prm;

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
  f_real *f_fields[NR_FIELDS];

  // C data structures
  void *c_ctx;

  // did we allocate the fields / particles (otherwise, Fortran did)
  bool allocated;
};

// we keep this info global for now.

extern struct psc psc;

void psc_create();
void psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part);
void psc_destroy();

void psc_setup_parameters();
void psc_setup_fields_zero();
void psc_setup_fields_1();
void psc_setup_particles_1();
void psc_setup_particles_random_yz();
void psc_dump_particles(const char *fname);
void psc_save_particles_ref();
void psc_save_fields_ref();
void psc_check_currents_ref();
void psc_check_particles_ref();
void psc_check_fields_ref();
void psc_create_test_1(const char *ops_name);
void psc_create_test_2(const char *ops_name);
void psc_create_test_3(const char *ops_name);

void psc_push_part_yz();
void psc_push_part_z();
void psc_push_part_yz_a();
void psc_push_part_yz_b();

// various implementations of the psc
// (something like Fortran, generic C, CUDA, ...)

extern struct psc_ops psc_ops_fortran;
extern struct psc_ops psc_ops_generic_c;
extern struct psc_ops psc_ops_sse2; //Intel SIMD instructions

// Wrappers for Fortran functions
void PIC_push_part_yz();
void PIC_push_part_z();
void PIC_push_part_yz_a();
void PIC_push_part_yz_b();

// ----------------------------------------------------------------------
// other bits and hacks...

#define sqr(a) ((a) * (a))

#define HERE printf("HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__)

#endif
