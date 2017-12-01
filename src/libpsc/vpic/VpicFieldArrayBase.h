
#ifndef VPIC_FIELD_ARRAY_BASE_H
#define VPIC_FIELD_ARRAY_BASE_H

#include "grid.h"
#include "material.h"

#include "field_advance/field_advance.h"

#include <mrc_common.h>
#include <cassert>

// ======================================================================
// VpicFieldArrayBase

inline void field_array_ctor(field_array_t *fa, grid_t *g, const material_t *m_list, float damp);
inline void field_array_dtor(field_array_t *fa);

template<class ML>
struct VpicFieldArrayBase : field_array_t {
  typedef ML MaterialList;
  typedef field_t Element;
  
  enum {
    EX  = 0,
    EY  = 1,
    EZ  = 2,
    CBX = 4,
    CBY = 5,
    CBZ = 6,
    N_COMP = sizeof(field_t) / sizeof(float),
  };
  
  VpicFieldArrayBase(Grid* grid, MaterialList material_list, float damp)
  {
    field_array_ctor(this, grid, material_list, damp);
  }
  
  ~VpicFieldArrayBase()
  {
    field_array_dtor(this);
  }

  Element* data()
  {
    return f;
  }
  
  float* getData(int* ib, int* im)
  {
    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    im[0] = g->nx + 2*B;
    im[1] = g->ny + 2*B;
    im[2] = g->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &f[0].ex;
  }

  // These operators can be used to access the field directly,
  // though the performance isn't great, so one you use Field3D
  // when performance is important
  float operator()(int m, int i, int j, int k) const
  {
    float *f_ = &f[0].ex;
    return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
  }
  
  float& operator()(int m, int i, int j, int k)
  {
    float *f_ = &f[0].ex;
    return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
  }

  Element operator[](int idx) const
  {
    return f[idx];
  }
  
  Element& operator[](int idx)
  {
    return f[idx];
  }

  Grid* getGrid()
  {
    return static_cast<Grid*>(g);
  }
  
  // I'm keeping these for now, because I tink they're a nice interface,
  // but it doesn't scale well to other kinds of fields (as one can tell
  // from the macro use...)
#define MK_COMP_ACCESSOR(cbx)			\
  float cbx(int i, int j, int k) const		\
  {						\
    const int nx = g->nx, ny = g->ny;		\
    return f[VOXEL(i,j,k, nx,ny,nz)].cbx;	\
  }						\
						\
  float& cbx(int i, int j, int k)		\
  {						\
    const int nx = g->nx, ny = g->ny;		\
    return f[VOXEL(i,j,k, nx,ny,nz)].cbx;	\
  }

  MK_COMP_ACCESSOR(cbx)
  MK_COMP_ACCESSOR(cby)
  MK_COMP_ACCESSOR(cbz)
  MK_COMP_ACCESSOR(ex)
  MK_COMP_ACCESSOR(ey)
  MK_COMP_ACCESSOR(ez)
  
};

// ----------------------------------------------------------------------
// copy&paste from sfa.c -- the only reason I can't just use this
// is that there's a new/delete interface, but no ctor/dtor interface

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

inline void
destroy_sfa_params( sfa_params_t * p ) {
  FREE_ALIGNED( p->mc );
  FREE( p );
}

// ----------------------------------------------------------------------
// the following are basically the same as the corresponding new/delete
// functions from sfa.c, but changed to leave out the allocation / free

inline void field_array_ctor(field_array_t *fa, grid_t *g, const material_t *m_list, float damp)
{
  assert(g && m_list && damp >= 0.);
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
}

inline void field_array_dtor(field_array_t *fa)
{
  destroy_sfa_params( (sfa_params_t *)fa->params );
  FREE_ALIGNED( fa->f );
}

// ======================================================================
// Field3D
//
// A class to accelerate 3-d field access
// It works for VpicFieldArray and VpicInterpolatorArray, but should
// be generalized to work without knowing their internals

template<class FA>
struct Field3D {
  typedef FA Array;
  typedef typename FA::Element Element;
  
  Field3D(Array& fa)
    : sx_(fa.g->nx + 2), sy_(fa.g->ny + 2),
      f_(fa.data())
  {
  }

  int voxel(int i, int j, int k) const
  {
    return i + sx_ * (j + sy_ * (k));
  }

  // FIXME, this is an odd mix of two interfaces
  // first field3d just being used for grid information,
  // and providing access to the whole struct
  // This interface can't easily be converted to SOA
  Element& operator()(Array &fa, int i, int j, int k)
  {
    return fa.data()[voxel(i,j,k)];
  }
  
  Element operator()(Array &fa, int i, int j, int k) const
  {
    return fa.data()[voxel(i,j,k)];
  }

  Element& operator()(int i, int j, int k)
  {
    return f_[voxel(i,j,k)];
  }

  Element operator()(int i, int j, int k) const
  {
    return f_[voxel(i,j,k)];
  }

  // second, access to a particular component, but this one is
  // for the specific Array used at construction time
  float& operator()(int m, int i, int j, int k)
  {
    float * RESTRICT f = reinterpret_cast<float *>(f_);
    return f[m + Array::N_COMP * voxel(i,j,k)];
  }
  
  float operator()(int m, int i, int j, int k) const
  {
    float * RESTRICT f = reinterpret_cast<float *>(f_);
    return f[m + Array::N_COMP * voxel(i,j,k)];
  }

private:
  int sx_, sy_;
  Element * RESTRICT f_;
};



#endif

