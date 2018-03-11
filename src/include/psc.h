/*! @file */

#ifndef PSC_H
#define PSC_H

#include "psc_config.h"

#include "psc_bits.h"
#include "psc_particles.h"
#include "psc_fields.h"
#include "grid.hxx"

#include <mrc_domain.h>

#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <vector>

// ----------------------------------------------------------------------
// cell_map

struct cell_map {
  int dims[3];
  int b_bits[3];
  int b_bits_max;
  int N; // indices will be 0..N-1
};

int cell_map_init(struct cell_map *map, const int dims[3],
		  const int blocksize[3]);
int cell_map_3to1(struct cell_map *map, int i[3]);
void cell_map_1to3(struct cell_map *map, int idx, int i[3]);
void cell_map_free(struct cell_map *map);

// ----------------------------------------------------------------------

enum {
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
  NR_FIELDS,
};

// C floating point type
// used to switch between single and double precision

typedef float real;

#define real(x) x ## f

// Fortran types

typedef double f_real;
typedef int f_int;

///User specified parameters
///
struct psc_param {
  double qq;	///<elemental charge 
  double mm;	///<mass
  double tt;	///<some measurement for energy ? (default is 1keV in fortran) 
  double cc;	///<speed of light
  double eps0;	///<vacuum permittivity
  int nmax;	///<number of timesteps
  double lw;	///<normalization coefficient for laser wavelength (omega)
  double i0;	///<laser intensity
  double n0;	///<electron density
  double e0;	///<field intesity
  double b0;
  double j0;
  double rho0;
  double phi0;
  double a0;
  double cfl;   ///<CFL number to be used for determining timestep
  int nicell;	///<number of particles per gridpoint to represent a normalised density of 1 
  int nr_populations;  ///< nr of different particle populations (defaults to nr_kinds)
  int neutralizing_population;  ///< the initial number of particles in a cell for this population will be st so that it achieves neutrality
  bool seed_by_time;
  bool fractional_n_particles_per_cell;
  bool const_num_particles_per_cell;
  bool initial_momentum_gamma_correction;
  double wallclock_limit;
  bool write_checkpoint;
  int write_checkpoint_every_step;
  char *fields_base; ///< base type for psc_mfields ("c", "fortran", "cuda")
  char *particles_base; ///< base type for psc_mparticles ("c", "fortran", "cuda")
  int stats_every; ///< output timing and other info every so many steps
  bool detailed_profiling; ///< output profiling info for each process separately
  double theta_xz; ///< rotate anisotropic maxwellian in x-z plane

  int sort_interval;
};

/// coefficients needed for computations
/// -- derived, not provided by user
struct psc_coeff {
  double cori;	///< 1 / psc_params.nicell 
  double alpha;
  double beta;
  double eta;

  // FIXME are these needed in general?
  double wl;	///<omega
  double ld;	///Normalization factor for lengths
  double vos;	///< 1/k
  double vt;
  double wp;
};

// need to match fortran values

///Possible boundary conditions for fields
enum {
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_UPML,
  BND_FLD_TIME,
  BND_FLD_CONDUCTING_WALL,
  BND_FLD_ABSORBING,
};

///Possible boundary conditions for particles
enum {
  BND_PART_REFLECTING,
  BND_PART_PERIODIC,
  BND_PART_ABSORBING,
  BND_PART_OPEN,
};

///Describes the spatial domain to operate on.
///
///This struct describes the spatial dimension of the simulation-box
///@note Here, you can also set the dimensionality by eliminating a dimension. Example: To simulate in xy only, set
///\verbatim psc_domain.gdims[2]=1 \endverbatim
///Also, set the boundary conditions for the eliminated dimensions to BND_FLD_PERIODIC or you'll get invalid \a dt and \a dx

struct psc_domain {
  double length[3];	///<The physical size of the simulation-box 
  double corner[3];
  int gdims[3];		///<Number of grid-points in each dimension
  int np[3];		///<Number of patches in each dimension
  int bs[3];
  int bnd_fld_lo[3];	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  int bnd_fld_hi[3];	///<Boundary conditions of the fields. Can be any value of BND_FLD.
  int bnd_part_lo[3];	///<Boundary conditions of the particles. Can be any value of BND_PART.
  int bnd_part_hi[3];   ///<Boundary conditions of the particles. Can be any value of BND_PART.
};

// ----------------------------------------------------------------------
// general info / parameters for the code

///Describes the different particle kinds
///
struct psc_kind {
  double q;  // charge
  double m;  // mass
  double n;  // default density
  double T;  // default temperature
  char *name; // short string for the kind name
};

///Default kinds (electrons + ions)
enum {
  KIND_ELECTRON,
  KIND_ION,
  NR_KINDS,
};

#define CRDX(p, jx) (psc->grid().dx[0] * (jx) + psc->grid().patches[p].xb[0])
#define CRDY(p, jy) (psc->grid().dx[1] * (jy) + psc->grid().patches[p].xb[1])
#define CRDZ(p, jz) (psc->grid().dx[2] * (jz) + psc->grid().patches[p].xb[2])

///This structure holds all the interfaces for the given configuration.
///
///
struct psc {
  ///@defgroup interfaces Interfaces @{
  struct mrc_obj obj;
  struct psc_method *method;                    ///< particular variant of PIC method
  struct psc_push_particles *push_particles;	///< particle pusher
  struct psc_push_fields *push_fields;		///< field pusher
  struct psc_bnd *bnd;				///< boundaries
  struct psc_bnd_particles *bnd_particles;	///< boundary particlesxs
  struct psc_collision *collision;		///< collision operator
  struct psc_sort *sort;			///< sort operator
  struct psc_marder *marder;                    ///< marder correction
  struct psc_diag *diag;                	///< timeseries diagnostics
  struct psc_output_fields_collection *output_fields_collection; ///< collection of psc_output_fields
  struct psc_output_particles *output_particles;///< particle output
  struct psc_event_generator *event_generator;	///< event generator
  struct psc_balance *balance;                  ///< rebalancer
  struct psc_checks *checks;                    ///< run-time checks
  ///@}

  ///@defgroup config-params user-configurable parameters @{
  // user-configurable parameters
  struct psc_param prm;		///< normalization parameters set by the user
  struct psc_coeff coeff;	///< automatically derived constants
  struct psc_domain domain;	///< the computational domain
  int nr_kinds;                 ///< nr of different particle kinds
  struct psc_kind *kinds;       ///< particle kinds (e.g., e-, ion, ...)
  ///@}

  
  // other parameters / constants
  double p2A, p2B;
  int timestep;	///< the current timestep
  double dt;	///< timestep in physical units
  ///@}

  struct psc_mparticles *particles;	///< All the particles, indexed by their containing patch
  struct psc_mfields *flds;	///< The fields.
  int n_state_fields;           ///< How many field components do we need in ::flds

  ///The domain partitioner.
  ///
  ///Use this to access the global list of patches \sa \ref patches

  struct mrc_domain *mrc_domain;

  Grid_t* make_grid(struct mrc_domain* domain);
  int n_patches() { return grid().n_patches(); }

  int ibn[3];         ///< number of ghost points

  double time_start;

  const Grid_t& grid() const { assert(grid_); return *grid_; }

  Grid_t* grid_ = {};
};

MRC_CLASS_DECLARE(psc, struct psc);

struct psc_particle_npt {
  int kind; ///< particle kind
  double q; ///< charge
  double m; ///< mass
  double n; ///< density
  double p[3]; ///< momentum
  double T[3]; ///< temperature
  int particles_per_cell; ///< desired number of particles per cell per unit density. If not specified, the global nicell is used.
};

struct psc_ops {
  MRC_SUBCLASS_OPS(struct psc);
  // FIXME setup_particles -> set_ic_particles
  void (*setup_particles)(struct psc *psc, uint *nr_particles_by_patch, bool count_only);
  void (*init_npt)(struct psc *psc, int kind, double x[3],
		   struct psc_particle_npt *npt);
  // FIXME setup_fields -> set_ic_fields
  void (*setup_fields)(struct psc *psc, struct psc_mfields *flds);
  double (*init_field)(struct psc *psc, double x[3], int m);
  void (*integrate)(struct psc *psc);
  void (*step)(struct psc *psc);
  void (*output)(struct psc *psc);
  struct mrc_domain *(*setup_mrc_domain)(struct psc *psc, int nr_patches);
};

#define psc_ops(psc) ((struct psc_ops *)((psc)->obj.ops))

/*!
Enumerates all grid points of patch p

Use this macro like a regular for() expression to process all grid points on any patch p.
Variables ix, iy, iz will be automatically declared.

\param p the patch to process
\param l Number of ghost-points to include in negative direction
\param r Number of ghost-points to include in positive direction
\param[out] ix x-coordinate of current grid point
\param[out] iy y-coordinate of current grid point
\param[out] iz z-coordinate of current grid point

Always close this expression with foreach_3d_end
*/
#define foreach_3d(psc, p, ix, iy, iz, l, r) {                         \
  int __ilo[3] = { (psc)->domain.gdims[0] == 1 ? 0 : -l ,              \
                  (psc)->domain.gdims[1] == 1 ? 0 : -l,                        \
                  (psc)->domain.gdims[2] == 1 ? 0 : -l };              \
  int __ihi[3] = { (psc)->grid().ldims[0] + ((psc)->domain.gdims[0] == 1 ? 0 : r), \
		   (psc)->grid().ldims[1] + ((psc)->domain.gdims[1] == 1 ? 0 : r), \
		   (psc)->grid().ldims[2] + ((psc)->domain.gdims[2] == 1 ? 0 : r) }; \
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

/*! Closes a loop started by foreach_3d
 * 
Usage:
\code
psc_foreach_patch(&psc, p){
psc_foreach_3d(p, jx, jy, jz, 0, 0) {
  //do stuff
} psc_foreach_3d_end;
}
\endcode
*/
#define foreach_3d_end				\
  } } }

#define psc_foreach_3d(psc, p, ix, iy, iz, l, r) {			\
  int __ilo[3] = { -l, -l, -l };					\
  int __ihi[3] = { psc->grid().ldims[0] + r,				\
		   psc->grid().ldims[1] + r,				\
		   psc->grid().ldims[2] + r };				\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_end				\
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

#define psc_foreach_3d_g(psc, p, ix, iy, iz) {				\
  int __ilo[3] = { -psc->ibn[0], -psc->ibn[1], -psc->ibn[2] };		\
  int __ihi[3] = { psc->grid().ldims[0] + psc->ibn[0],			\
		   psc->grid().ldims[1] + psc->ibn[1],			\
		   psc->grid().ldims[2] + psc->ibn[2] };		\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {			\
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {			\
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_g_end				\
  } } }

#define psc_foreach_patch(psc, p)		\
  for (int p = 0; p < (psc)->n_patches(); p++)


// ----------------------------------------------------------------------
// we keep this info global for now.

extern struct psc *ppsc;

struct psc *psc_create(MPI_Comm comm);
void psc_set_from_options(struct psc *psc);
void psc_setup(struct psc *psc);
void psc_set_kinds(struct psc *psc, int nr_kinds, const struct psc_kind *kinds);
void psc_view(struct psc *psc);
void psc_destroy(struct psc *psc);
struct psc_particle_double;
void psc_setup_particle(struct psc *psc, struct psc_particle_double *prt, struct psc_particle_npt *npt,
			int p, double xx[3]);
void psc_setup_partition(struct psc *psc, uint *nr_particles_by_patch);
void psc_setup_particles(struct psc *psc, uint *nr_particles_by_patch);
void psc_set_ic_fields(struct psc *psc);
void psc_set_ic_fields_default(struct psc *psc);
void psc_output(struct psc *psc);
void psc_integrate(struct psc *psc);

void psc_setup_coeff(struct psc *psc);
void psc_setup_domain(struct psc *psc);
struct mrc_domain *psc_setup_mrc_domain(struct psc *psc, int nr_patches);

struct psc *psc_read_checkpoint(MPI_Comm comm, int n);
void psc_write_checkpoint(struct psc *psc);

void psc_setup_fortran(struct psc *psc);
void psc_print_profiling(struct psc *psc);


void psc_default_dimensionless(struct psc *psc);

int psc_main(int *argc, char ***argv, struct psc_ops *type);

static inline bool psc_at_boundary_lo(struct psc *psc, int p, int d)
{
  return psc->grid().patches[p].off[d] == 0;
}

static inline bool psc_at_boundary_hi(struct psc *psc, int p, int d)
{
  return psc->grid().patches[p].off[d] + psc->grid().ldims[d] == psc->grid().gdims[d];
}

// ----------------------------------------------------------------------
// psc_stats: simple statistics

#define MAX_PSC_STATS 20

extern double psc_stats_val[MAX_PSC_STATS+1]; // [0] is left empty
extern int nr_psc_stats;

int psc_stats_register(const char *name);
void psc_stats_log(struct psc *psc);

#define psc_stats_start(n) do {				\
    psc_stats_val[n] -= MPI_Wtime();			\
  } while (0)

#define psc_stats_stop(n) do {				\
    psc_stats_val[n] += MPI_Wtime();			\
  } while (0)

// These are general statistics categories to be used in different parts
// of the code as appropriate.

extern int st_time_output;   //< time spent in output
extern int st_time_comm;     //< time spent in communications
extern int st_time_particle; //< time spent in particle computation
extern int st_time_field;    //< time spent in field computation

// ----------------------------------------------------------------------
// other bits and hacks...

#endif
