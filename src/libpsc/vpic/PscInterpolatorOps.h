
#ifndef PSC_INTERPOLATOR_OPS_H
#define PSC_INTERPOLATOR_OPS_H

#define fi(x,y,z) fi[   VOXEL(x,y,z, nx,ny,nz) ]
#define f(x,y,z)  f [   VOXEL(x,y,z, nx,ny,nz) ]
#define nb(x,y,z) nb[ 6*VOXEL(x,y,z, nx,ny,nz) ]

typedef struct load_interpolator_pipeline_args {
  MEM_PTR( interpolator_t, 128 ) fi;
  MEM_PTR( const field_t,  128 ) f;
  MEM_PTR( const int64_t,  128 ) nb;
  int nx;
  int ny;
  int nz;
} load_interpolator_pipeline_args_t;

static inline void
load_interpolator(interpolator_array_t* ia, const field_array_t* fa)
{
  interpolator_t * ALIGNED(128) fi = ia->i;
  const field_t  * ALIGNED(128) f  = fa->f;

  interpolator_t * ALIGNED(16) pi;

  const field_t  * ALIGNED(16) pf0;
  const field_t  * ALIGNED(16) pfx,  * ALIGNED(16) pfy,  * ALIGNED(16) pfz;
  const field_t  * ALIGNED(16) pfyz, * ALIGNED(16) pfzx, * ALIGNED(16) pfxy;
  int x, y, z, n_voxel;

  const int nx = ia->g->nx;
  const int ny = ia->g->ny;
  const int nz = ia->g->nz;

  const float fourth = 0.25;
  const float half   = 0.5;

  float w0, w1, w2, w3;

  // Process the voxels assigned to this pipeline
  
  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 1,
                     0, 1, x, y, z, n_voxel );

# define LOAD_STENCIL()    \
  pi   = &fi(x,  y,  z  ); \
  pf0  =  &f(x,  y,  z  ); \
  pfx  =  &f(x+1,y,  z  ); \
  pfy  =  &f(x,  y+1,z  ); \
  pfz  =  &f(x,  y,  z+1); \
  pfyz =  &f(x,  y+1,z+1); \
  pfzx =  &f(x+1,y,  z+1); \
  pfxy =  &f(x+1,y+1,z  )

  for( ; n_voxel; n_voxel-- ) {
    LOAD_STENCIL();

    // ex interpolation
    w0 = pf0->ex;
    w1 = pfy->ex;
    w2 = pfz->ex;
    w3 = pfyz->ex;
    pi->ex       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dexdy    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dexdz    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2exdydz = fourth*( (w3 + w0) - (w1 + w2) );

    // ey interpolation coefficients
    w0 = pf0->ey;
    w1 = pfz->ey;
    w2 = pfx->ey;
    w3 = pfzx->ey;
    pi->ey       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->deydz    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->deydx    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2eydzdx = fourth*( (w3 + w0) - (w1 + w2) );

    // ez interpolation coefficients
    w0 = pf0->ez;
    w1 = pfx->ez;
    w2 = pfy->ez;
    w3 = pfxy->ez;
    pi->ez       = fourth*( (w3 + w0) + (w1 + w2) );
    pi->dezdx    = fourth*( (w3 - w0) + (w1 - w2) );
    pi->dezdy    = fourth*( (w3 - w0) - (w1 - w2) );
    pi->d2ezdxdy = fourth*( (w3 + w0) - (w1 + w2) );

    // bx interpolation coefficients
    w0 = pf0->cbx;
    w1 = pfx->cbx;
    pi->cbx    = half*( w1 + w0 );
    pi->dcbxdx = half*( w1 - w0 );

    // by interpolation coefficients
    w0 = pf0->cby;
    w1 = pfy->cby;
    pi->cby    = half*( w1 + w0 );
    pi->dcbydy = half*( w1 - w0 );

    // bz interpolation coefficients
    w0 = pf0->cbz;
    w1 = pfz->cbz;
    pi->cbz    = half*( w1 + w0 );
    pi->dcbzdz = half*( w1 - w0 );

    pi++; pf0++; pfx++; pfy++; pfz++; pfyz++; pfzx++; pfxy++;

    x++;
    if( x>nx ) {
      x=1, y++;
      if( y>ny ) y=1, z++;
    }
  }

# undef LOAD_STENCIL

}

template<class I, class F>
struct PscInterpolatorOps {
  typedef I Interpolator;
  typedef F FieldArray;

  void load_interpolator_array(Interpolator *interpolator,
			       FieldArray *vmflds)
  {
    TIC ::load_interpolator(interpolator, vmflds); TOC(load_interpolator, 1);
  }
};

#endif

