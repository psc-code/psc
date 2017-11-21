
#ifndef PSC_INTERPOLATOR_OPS_H
#define PSC_INTERPOLATOR_OPS_H

#define fi(x,y,z) fi[   VOXEL(x,y,z, nx,ny,nz) ]
#define f(x,y,z)  f [   VOXEL(x,y,z, nx,ny,nz) ]

static inline void
load_interpolator(interpolator_array_t* ia, const field_array_t* fa)
{
  interpolator_t * ALIGNED(128) fi = ia->i;
  const field_t  * ALIGNED(128) f  = fa->f;

  const int nx = ia->g->nx;
  const int ny = ia->g->ny;
  const int nz = ia->g->nz;

  const float fourth = 0.25;
  const float half   = 0.5;

  float w0, w1, w2, w3;

  for (int z = 1; z <= nz; z++) {
    for (int y = 1; y <= ny; y++) {
      for (int x = 1; x <= nx; x++) {

    // ex interpolation
    w0 = f(x,y  ,z  ).ex;
    w1 = f(x,y+1,z  ).ex;
    w2 = f(x,y  ,z+1).ex;
    w3 = f(x,y+1,z+1).ex;
    fi(x,y,z).ex       = fourth*( (w3 + w0) + (w1 + w2) );
    fi(x,y,z).dexdy    = fourth*( (w3 - w0) + (w1 - w2) );
    fi(x,y,z).dexdz    = fourth*( (w3 - w0) - (w1 - w2) );
    fi(x,y,z).d2exdydz = fourth*( (w3 + w0) - (w1 + w2) );

    // ey interpolation coefficients
    w0 = f(x  ,y,z  ).ey;
    w1 = f(x  ,y,z+1).ey;
    w2 = f(x+1,y,z  ).ey;
    w3 = f(x+1,y,z+1).ey;
    fi(x,y,z).ey       = fourth*( (w3 + w0) + (w1 + w2) );
    fi(x,y,z).deydz    = fourth*( (w3 - w0) + (w1 - w2) );
    fi(x,y,z).deydx    = fourth*( (w3 - w0) - (w1 - w2) );
    fi(x,y,z).d2eydzdx = fourth*( (w3 + w0) - (w1 + w2) );

    // ez interpolation coefficients
    w0 = f(x  ,y,z  ).ez;
    w1 = f(x+1,y,z  ).ez;
    w2 = f(x  ,y+1,z).ez;
    w3 = f(x+1,y+1,z).ez;
    fi(x,y,z).ez       = fourth*( (w3 + w0) + (w1 + w2) );
    fi(x,y,z).dezdx    = fourth*( (w3 - w0) + (w1 - w2) );
    fi(x,y,z).dezdy    = fourth*( (w3 - w0) - (w1 - w2) );
    fi(x,y,z).d2ezdxdy = fourth*( (w3 + w0) - (w1 + w2) );

    // bx interpolation coefficients
    w0 = f(x  ,y,z).cbx;
    w1 = f(x+1,y,z).cbx;
    fi(x,y,z).cbx    = half*( w1 + w0 );
    fi(x,y,z).dcbxdx = half*( w1 - w0 );

    // by interpolation coefficients
    w0 = f(x,y  ,z).cby;
    w1 = f(x,y+1,z).cby;
    fi(x,y,z).cby    = half*( w1 + w0 );
    fi(x,y,z).dcbydy = half*( w1 - w0 );

    // bz interpolation coefficients
    w0 = f(x,y,z  ).cbz;
    w1 = f(x,y,z+1).cbz;
    fi(x,y,z).cbz    = half*( w1 + w0 );
    fi(x,y,z).dcbzdz = half*( w1 - w0 );
      }
    }
  }
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

