
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
  void mp_size_recv_buffer(int tag, int size);
  void mp_size_send_buffer(int tag, int size);

  //private:
  grid *g_;
};

inline Grid::Grid(grid *g)
  : g_(g)
{
}

inline void Grid::setup(double dx[3], double dt, double cvac, double eps0)
{
  g_->dx = dx[0];
  g_->dy = dx[1];
  g_->dz = dx[2];
  g_->dt = dt;
  g_->cvac = cvac;
  g_->eps0 = eps0;
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

inline void Grid::mp_size_recv_buffer(int tag, int size)
{
  ::mp_size_recv_buffer(g_->mp, tag, size);
}

inline void Grid::mp_size_send_buffer(int tag, int size)
{
  ::mp_size_send_buffer(g_->mp, tag, size);
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

struct FieldArray : field_array_t {
};

// inline FieldArrayPtr::FieldArrayPtr(field_array_t*& p)
//   : p_(p)
// {
// }

// ======================================================================
// FieldArrayPtr

struct FieldArrayPtr {
  FieldArrayPtr(field_array_t*& p);

  field_array_t *&p_;
};

inline FieldArrayPtr::FieldArrayPtr(field_array_t*& p)
  : p_(p)
{
}

// ======================================================================
// InterpolatorArrayPtr

struct InterpolatorArrayPtr {
  InterpolatorArrayPtr(interpolator_array_t*& p);

  interpolator_array_t *&p_;
};

inline InterpolatorArrayPtr::InterpolatorArrayPtr(interpolator_array_t*& p)
  : p_(p)
{
}

// ======================================================================
// AccumulatorArrayPtr

struct AccumulatorArrayPtr {
  AccumulatorArrayPtr(accumulator_array_t*& p);

  accumulator_array_t *&p_;
};

inline AccumulatorArrayPtr::AccumulatorArrayPtr(accumulator_array_t*& p)
  : p_(p)
{
}

// ======================================================================
// HydroArrayPtr

struct HydroArrayPtr {
  HydroArrayPtr(hydro_array_t*& p);

  hydro_array_t *&p_;
};

inline HydroArrayPtr::HydroArrayPtr(hydro_array_t*& p)
  : p_(p)
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
  field_array_t *new_field_array(double damp=0.);
  void define_field_array(double damp);
  
  RngPool rng_pool;

  //private:
  globals_diag *pDiag_;
  Grid grid_;
  MaterialList material_list_;
  FieldArrayPtr field_array_ptr_;
  InterpolatorArrayPtr interpolator_array_ptr_;
  AccumulatorArrayPtr accumulator_array_ptr_;
  HydroArrayPtr hydro_array_ptr_;
};

inline Simulation::Simulation()
  : grid_(simulation->grid),
    material_list_(simulation->material_list),
    field_array_ptr_(simulation->field_array),
    interpolator_array_ptr_(simulation->interpolator_array),
    accumulator_array_ptr_(simulation->accumulator_array),
    hydro_array_ptr_(simulation->hydro_array)
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

#define IN_sfa
#include "field_advance/standard/sfa_private.h"

static field_advance_kernels_t sfa_kernels = {

  // Destructor

  delete_standard_field_array,

  // Time stepping interfaces

  advance_b,
  advance_e,

  // Diagnostic interfaces

  energy_f,

  // Accumulator interfaces

  clear_jf,   synchronize_jf,
  clear_rhof, synchronize_rho,

  // Initialize interface

  compute_rhob,
  compute_curl_b,

  // Shared face cleaning interface

  synchronize_tang_e_norm_b,
  
  // Electric field divergence cleaning interface

  compute_div_e_err,
  compute_rms_div_e_err,
  clean_div_e,

  // Magnetic field divergence cleaning interface

  compute_div_b_err,
  compute_rms_div_b_err,
  clean_div_b

};

static float
minf( float a, 
      float b ) {
  return a<b ? a : b;
}

static sfa_params_t *
create_sfa_params( grid_t           * g,
                   const material_t * m_list,
                   float              damp ) {
  sfa_params_t * p;
  float ax, ay, az, cg2;
  material_coefficient_t *mc;
  const material_t *m;
  int n_mc;

  // Run sanity checks on the material list

  ax = g->nx>1 ? g->cvac*g->dt*g->rdx : 0; ax *= ax;
  ay = g->ny>1 ? g->cvac*g->dt*g->rdy : 0; ay *= ay;
  az = g->nz>1 ? g->cvac*g->dt*g->rdz : 0; az *= az;
  n_mc = 0;
  LIST_FOR_EACH(m,m_list) {
    if( m->sigmax/m->epsx<0 )
      WARNING(("\"%s\" is an active medium along x", m->name));
    if( m->epsy*m->muz<0 )
      WARNING(("\"%s\" has an imaginary x speed of light (ey)", m->name));
    if( m->epsz*m->muy<0 )
      WARNING(("\"%s\" has an imaginary x speed of light (ez)", m->name));
    if( m->sigmay/m->epsy<0 )
      WARNING(("\"%s\" is an active medium along y", m->name));
    if( m->epsz*m->mux<0 )
      WARNING(("\"%s\" has an imaginary y speed of light (ez)", m->name));
    if( m->epsx*m->muz<0 )
      WARNING(("\"%s\" has an imaginary y speed of light (ex)", m->name));
    if( m->sigmaz/m->epsz<0 )
      WARNING(("\"%s\" is an an active medium along z", m->name));
    if( m->epsx*m->muy<0 )
      WARNING(("\"%s\" has an imaginary z speed of light (ex)", m->name));
    if( m->epsy*m->mux<0 )
      WARNING(("\"%s\" has an imaginary z speed of light (ey)", m->name));
    cg2 = ax/minf(m->epsy*m->muz,m->epsz*m->muy) +
          ay/minf(m->epsz*m->mux,m->epsx*m->muz) +
          az/minf(m->epsx*m->muy,m->epsy*m->mux);
    if( cg2>=1 )
      WARNING(( "\"%s\" Courant condition estimate = %e", m->name, sqrt(cg2) ));
    if( m->zetax!=0 || m->zetay!=0 || m->zetaz!=0 )
      WARNING(( "\"%s\" magnetic conductivity is not supported" ));
    n_mc++;
  }

  // Allocate the sfa parameters

  MALLOC( p, 1 );
  MALLOC_ALIGNED( p->mc, n_mc+2, 128 );
  p->n_mc = n_mc;
  p->damp = damp;

  // Fill up the material coefficient array
  // FIXME: THIS IMPLICITLY ASSUMES MATERIALS ARE NUMBERED CONSECUTIVELY FROM
  // O.

  LIST_FOR_EACH( m, m_list ) {
    mc = p->mc + m->id;

    // Advance E coefficients
    // Note: m ->sigma{x,y,z} = 0 -> Non conductive
    //       mc->decay{x,y,z} = 0 -> Perfect conductor to numerical precision
    //       otherwise            -> Conductive
    ax = ( m->sigmax*g->dt ) / ( m->epsx*g->eps0 );
    ay = ( m->sigmay*g->dt ) / ( m->epsy*g->eps0 );
    az = ( m->sigmaz*g->dt ) / ( m->epsz*g->eps0 );
    mc->decayx = exp(-ax);
    mc->decayy = exp(-ay);
    mc->decayz = exp(-az);
    if( ax==0 )              mc->drivex = 1./m->epsx;
    else if( mc->decayx==0 ) mc->drivex = 0;
    else mc->drivex = 2.*exp(-0.5*ax)*sinh(0.5*ax) / (ax*m->epsx);
    if( ay==0 )              mc->drivey = 1./m->epsy;
    else if( mc->decayy==0 ) mc->drivey = 0;
    else mc->drivey = 2.*exp(-0.5*ay)*sinh(0.5*ay) / (ay*m->epsy);
    if( az==0 )              mc->drivez = 1./m->epsz;
    else if( mc->decayz==0 ) mc->drivez = 0;
    else mc->drivez = 2.*exp(-0.5*az)*sinh(0.5*az) / (az*m->epsz);
    mc->rmux = 1./m->mux;
    mc->rmuy = 1./m->muy;
    mc->rmuz = 1./m->muz;

    // Clean div E coefficients.  Note: The charge density due to J =
    // sigma E currents is not computed.  Consequently, the divergence
    // error inside conductors cannot computed.  The divergence error
    // multiplier is thus set to zero to ignore divergence errors
    // inside conducting materials.

    mc->nonconductive = ( ax==0 && ay==0 && az==0 ) ? 1. : 0.;
    mc->epsx = m->epsx;
    mc->epsy = m->epsy;
    mc->epsz = m->epsz;
  }

  return p;
}

inline field_array_t *Simulation::new_field_array(double damp)
{
  //return ::new_standard_field_array(grid_.g_, material_list_.ml_, damp);
  grid_t *g = grid_.g_;
  material_t *m_list = material_list_.ml_;
  
  field_array_t * fa;
  if( !g || !m_list || damp<0 ) ERROR(( "Bad args" ));
  MALLOC( fa, 1 );
  MALLOC_ALIGNED( fa->f, g->nv, 128 );
  CLEAR( fa->f, g->nv );
  fa->g = g;
  fa->params = create_sfa_params( g, m_list, damp );
  fa->kernel[0] = sfa_kernels;
  if( !m_list->next ) {
    /* If there is only one material, then this material permeates all
       space and we can use high performance versions of some kernels. */
    fa->kernel->advance_e         = vacuum_advance_e;
    fa->kernel->energy_f          = vacuum_energy_f;
    fa->kernel->compute_rhob      = vacuum_compute_rhob;
    fa->kernel->compute_curl_b    = vacuum_compute_curl_b;
    fa->kernel->compute_div_e_err = vacuum_compute_div_e_err;
    fa->kernel->clean_div_e       = vacuum_clean_div_e;
  }

  return fa;
}

inline void Simulation::define_field_array(double damp)
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
  
  field_array_ptr_.p_ = new_field_array(damp);
  interpolator_array_ptr_.p_ = ::new_interpolator_array(grid);
  accumulator_array_ptr_.p_ = ::new_accumulator_array(grid);
  hydro_array_ptr_.p_ = ::new_hydro_array(grid);
 
  // Pre-size communications buffers. This is done to get most memory
  // allocation over with before the simulation starts running
  
  int nx1 = grid->nx+1, ny1 = grid->ny+1, nz1 = grid->nz+1;
  grid_.mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
  grid_.mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
  grid_.mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
  grid_.mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
  grid_.mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
  grid_.mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
  
  grid_.mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(hydro_t));
  grid_.mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(hydro_t));
  grid_.mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(hydro_t));
  grid_.mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(hydro_t));
  grid_.mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(hydro_t));
  grid_.mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(hydro_t));
}


#endif



