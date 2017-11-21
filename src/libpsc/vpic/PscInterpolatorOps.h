
#ifndef PSC_INTERPOLATOR_OPS_H
#define PSC_INTERPOLATOR_OPS_H

#define fi(x,y,z) fi[   VOXEL(x,y,z, nx,ny,nz) ]
#define f(x,y,z)  f [   VOXEL(x,y,z, nx,ny,nz) ]

template<class I, class F>
struct PscInterpolatorOps {
  typedef I Interpolator;
  typedef F FieldArray;

  void load_interpolator(Interpolator* ia, const FieldArray* fa)
  {
    interpolator_t * ALIGNED(128) fi = ia->i;
    const field_t  * ALIGNED(128) f  = fa->f;
    
    const int nx = ia->g->nx;
    const int ny = ia->g->ny;
    const int nz = ia->g->nz;
    
    const float fourth = 0.25;
    const float half   = 0.5;
    
    float w0, w1, w2, w3;
    
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  
	  // ex interpolation
	  w0 = f(i,j  ,k  ).ex;
	  w1 = f(i,j+1,k  ).ex;
	  w2 = f(i,j  ,k+1).ex;
	  w3 = f(i,j+1,k+1).ex;
	  fi(i,j,k).ex       = fourth*((w3 + w0) + (w1 + w2));
	  fi(i,j,k).dexdy    = fourth*((w3 - w0) + (w1 - w2));
	  fi(i,j,k).dexdz    = fourth*((w3 - w0) - (w1 - w2));
	  fi(i,j,k).d2exdydz = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ey interpolation coefficients
	  w0 = f(i  ,j,k  ).ey;
	  w1 = f(i  ,j,k+1).ey;
	  w2 = f(i+1,j,k  ).ey;
	  w3 = f(i+1,j,k+1).ey;
	  fi(i,j,k).ey       = fourth*((w3 + w0) + (w1 + w2));
	  fi(i,j,k).deydz    = fourth*((w3 - w0) + (w1 - w2));
	  fi(i,j,k).deydx    = fourth*((w3 - w0) - (w1 - w2));
	  fi(i,j,k).d2eydzdx = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ez interpolation coefficients
	  w0 = f(i  ,j  ,k).ez;
	  w1 = f(i+1,j  ,k).ez;
	  w2 = f(i  ,j+1,k).ez;
	  w3 = f(i+1,j+1,k).ez;
	  fi(i,j,k).ez       = fourth*((w3 + w0) + (w1 + w2));
	  fi(i,j,k).dezdx    = fourth*((w3 - w0) + (w1 - w2));
	  fi(i,j,k).dezdy    = fourth*((w3 - w0) - (w1 - w2));
	  fi(i,j,k).d2ezdxdy = fourth*((w3 + w0) - (w1 + w2));
	  
	  // bx interpolation coefficients
	  w0 = f(i  ,j,k).cbx;
	  w1 = f(i+1,j,k).cbx;
	  fi(i,j,k).cbx    = half*(w1 + w0);
	  fi(i,j,k).dcbxdx = half*(w1 - w0);
	  
	  // by interpolation coefficients
	  w0 = f(i,j  ,k).cby;
	  w1 = f(i,j+1,k).cby;
	  fi(i,j,k).cby    = half*(w1 + w0);
	  fi(i,j,k).dcbydy = half*(w1 - w0);
	  
	  // bz interpolation coefficients
	  w0 = f(i,j,k  ).cbz;
	  w1 = f(i,j,k+1).cbz;
	  fi(i,j,k).cbz    = half*(w1 + w0);
	  fi(i,j,k).dcbzdz = half*(w1 - w0);
	}
      }
    }
  }
  
  void load_interpolator_array(Interpolator *interpolator,
			       FieldArray *vmflds)
  {
    TIC load_interpolator(interpolator, vmflds); TOC(load_interpolator, 1);
  }
};

#endif

