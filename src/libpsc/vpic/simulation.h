
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

#include <psc.h> // FIXME, only need the BND_* constants

// ======================================================================
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct Rng : rng {
  double drand();
  double drandn();
  double uniform(double lo, double hi);
  double normal(double mu, double sigma);
};

// ----------------------------------------------------------------------
// Rng implementation

inline double Rng::drand()
{
  return ::drand(this);
}

inline double Rng::drandn()
{
  return ::drandn(this);
}

inline double Rng::uniform(double lo, double hi)
{
  double dx = drand();
  return lo*(1.-dx) + hi*dx;
}

inline double Rng::normal(double mu, double sigma)
{
  return mu + sigma*drandn();
}

// ======================================================================
// class RngPool

struct RngPool {
  RngPool();
  
  void seed(int base, int which);
  Rng *operator[](int n);

  rng_pool *rng_pool_;
};

// ----------------------------------------------------------------------
// RngPool implementation

inline RngPool::RngPool()
{
  //rng_pool_ = simulation->entropy;
  int new_rng = 2;
  rng_pool_ = new_rng_pool(new_rng, 0, 0);
}

inline void RngPool::seed(int base, int which)
{
  seed_rng_pool(rng_pool_, base, which);
}

inline Rng* RngPool::operator[](int n)
{
  return reinterpret_cast<Rng *>(rng_pool_->rng[n]);
}

// ======================================================================
// Grid

struct Grid {
  Grid(grid *grid);
  
  void setup(double dx[3], double dt, double cvac, double eps0);
  void partition_periodic_box(double xl[3], double xh[3], int gdims[3], int np[3]);
  void set_fbc(int boundary, int fbc);
  void set_pbc(int boundary, int pbc);

  //private:
  grid *g_;
};

inline void Grid::setup(double dx[3], double dt, double cvac, double eps0)
{
  g_->dx = dx[0];
  g_->dy = dx[1];
  g_->dz = dx[2];
  g_->dt = dt;
  g_->cvac = cvac;
  g_->eps0 = eps0;
}

inline Grid::Grid(grid *g)
  : g_(g)
{
}

inline void Grid::partition_periodic_box(double xl[3], double xh[3],
					 int gdims[3], int np[3])
{
  ::partition_periodic_box(g_, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			   gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

inline void Grid::set_fbc(int boundary, int fbc)
{
  ::set_fbc(g_, boundary, fbc);
}

inline void Grid::set_pbc(int boundary, int pbc)
{
  ::set_pbc(g_, boundary, pbc);
}

// ======================================================================
// MaterialList

struct MaterialList {
  MaterialList(material_t*& m);

  material_t* append(material_t* m);
  bool empty();
  
  //private:
  material_t *&ml_;
};

inline MaterialList::MaterialList(material_t*& m)
  : ml_(m)
{
}

inline material_t* MaterialList::append(material_t* m)
{
  return ::append_material(m, &ml_);
}

inline bool MaterialList::empty()
{
  return !ml_;
}

// ======================================================================
// FieldArray

struct FieldArray {
  FieldArray(field_array_t*& a);

  field_array_t *&a_;
};

inline FieldArray::FieldArray(field_array_t*& a)
  : a_(a)
{
}

// ======================================================================
// InterpolatorArray

struct InterpolatorArray {
  InterpolatorArray(interpolator_array_t*& a);

  interpolator_array_t *&a_;
};

inline InterpolatorArray::InterpolatorArray(interpolator_array_t*& a)
  : a_(a)
{
}

// ======================================================================
// AccumulatorArray

struct AccumulatorArray {
  AccumulatorArray(accumulator_array_t*& a);

  accumulator_array_t *&a_;
};

inline AccumulatorArray::AccumulatorArray(accumulator_array_t*& a)
  : a_(a)
{
}

// ======================================================================
// HydroArray

struct HydroArray {
  HydroArray(hydro_array_t*& a);

  hydro_array_t *&a_;
};

inline HydroArray::HydroArray(hydro_array_t*& a)
  : a_(a)
{
}

// ======================================================================
// class Simulation

struct Simulation {
  Simulation();
  ~Simulation();

  void setup_grid(double dx[3], double dt, double cvac, double eps0);
  void define_periodic_grid(double xl[3], double xh[3], int gdims[3],
			    int np[3]);
  void set_domain_field_bc(int boundary, int bc);
  void set_domain_particle_bc(int boundary, int bc);

  struct material *define_material(const char *name, double eps, double mu,
				   double sigma, double zeta);
  inline void define_field_array(struct field_array *fa, double damp);
  
  RngPool rng_pool;

  //private:
  globals_diag *pDiag_;
  Grid grid_;
  MaterialList material_list_;
  FieldArray field_array_;
  InterpolatorArray interpolator_array_;
  AccumulatorArray accumulator_array_;
  HydroArray hydro_array_;
};

// ----------------------------------------------------------------------
// Simulation implementation

inline Simulation::Simulation()
  : grid_(simulation->grid),
    material_list_(simulation->material_list),
    field_array_(simulation->field_array),
    interpolator_array_(simulation->interpolator_array),
    accumulator_array_(simulation->accumulator_array),
    hydro_array_(simulation->hydro_array)
{
}

inline Simulation::~Simulation()
{
  delete pDiag_;
}

inline void Simulation::setup_grid(double dx[3], double dt, double cvac, double eps0)
{
  grid_.setup(dx, dt, cvac, eps0);
}

inline void Simulation::define_periodic_grid(double xl[3], double xh[3], int gdims[3],
					     int np[3])
{
  simulation->px = size_t(np[0]);
  simulation->py = size_t(np[1]);
  simulation->pz = size_t(np[2]);
  grid_.partition_periodic_box(xl, xh, gdims, np);
}

inline void Simulation::set_domain_field_bc(int boundary, int bc)
{
  int fbc;
  switch (bc) {
  case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
  case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
  default: assert(0);
  }
  grid_.set_fbc(boundary, fbc);
}

inline void Simulation::set_domain_particle_bc(int boundary, int bc)
{
  int pbc;
  switch (bc) {
  case BND_PART_REFLECTING: pbc = reflect_particles; break;
  case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
  default: assert(0);
  }
  grid_.set_pbc(boundary, pbc);
}

inline struct material *Simulation::define_material(const char *name,
						    double eps, double mu=1.,
						    double sigma=0., double zeta=0.)
{
  return material_list_.append(material(name,
				       eps,   eps,   eps,
				       mu,    mu,    mu,
				       sigma, sigma, sigma,
				       zeta,  zeta,  zeta));
}

inline void Simulation::define_field_array(struct field_array *fa, double damp)
{
  grid_t *grid = grid_.g_;
 
  if (grid->nx<1 || grid->ny<1 || grid->nz<1 ) {
    mprintf("Define your grid before defining the field array\n");
    assert(0);
  }
  if (material_list_.empty()) {
    mprintf("Define your materials before defining the field array\n");
    assert(0);
  }
  
  field_array_.a_ = fa ? fa :
    ::new_standard_field_array(grid, simulation->material_list, damp);
  interpolator_array_.a_ = new_interpolator_array(grid);
  accumulator_array_.a_ = new_accumulator_array(grid);
  hydro_array_.a_ = new_hydro_array(grid);
 
  // Pre-size communications buffers. This is done to get most memory
  // allocation over with before the simulation starts running
  
  int nx1 = grid->nx+1, ny1 = grid->ny+1, nz1 = grid->nz+1;
  mp_size_recv_buffer(grid->mp, BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
  mp_size_recv_buffer(grid->mp, BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
  mp_size_recv_buffer(grid->mp, BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
  mp_size_recv_buffer(grid->mp, BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
  mp_size_recv_buffer(grid->mp, BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
  mp_size_recv_buffer(grid->mp, BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  
  mp_size_send_buffer(grid->mp, BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
  mp_size_send_buffer(grid->mp, BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
  mp_size_send_buffer(grid->mp, BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
  mp_size_send_buffer(grid->mp, BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
  mp_size_send_buffer(grid->mp, BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
  mp_size_send_buffer(grid->mp, BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
}


#endif



