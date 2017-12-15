
#ifndef PSC_INTERPOLATOR_H
#define PSC_INTERPOLATOR_H

template<class InterpolatorBase, class FA>
struct PscInterpolator : InterpolatorBase
{
  typedef InterpolatorBase Base;
  typedef FA FieldArray;
  typedef typename Base::Grid Grid;

  using Base::g;

  static PscInterpolator* create(Grid *grid)
  {
    return reinterpret_cast<PscInterpolator*>(Base::create(grid));
  }
  
  // ----------------------------------------------------------------------
  // load
  
  void do_load(/*const*/ FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Field3D<Base> I(*this);

    const int nx = g->nx, ny = g->ny, nz = g->nz;
    const float fourth = 0.25;
    const float half   = 0.5;

    float w0, w1, w2, w3;
    
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  
	  // ex interpolation
	  w0 = F(i,j  ,k  ).ex;
	  w1 = F(i,j+1,k  ).ex;
	  w2 = F(i,j  ,k+1).ex;
	  w3 = F(i,j+1,k+1).ex;
	  I(i,j,k).ex       = fourth*((w3 + w0) + (w1 + w2));
	  I(i,j,k).dexdy    = fourth*((w3 - w0) + (w1 - w2));
	  I(i,j,k).dexdz    = fourth*((w3 - w0) - (w1 - w2));
	  I(i,j,k).d2exdydz = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ey interpolation coefficients
	  w0 = F(i  ,j,k  ).ey;
	  w1 = F(i  ,j,k+1).ey;
	  w2 = F(i+1,j,k  ).ey;
	  w3 = F(i+1,j,k+1).ey;
	  I(i,j,k).ey       = fourth*((w3 + w0) + (w1 + w2));
	  I(i,j,k).deydz    = fourth*((w3 - w0) + (w1 - w2));
	  I(i,j,k).deydx    = fourth*((w3 - w0) - (w1 - w2));
	  I(i,j,k).d2eydzdx = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ez interpolation coefficients
	  w0 = F(i  ,j  ,k).ez;
	  w1 = F(i+1,j  ,k).ez;
	  w2 = F(i  ,j+1,k).ez;
	  w3 = F(i+1,j+1,k).ez;
	  I(i,j,k).ez       = fourth*((w3 + w0) + (w1 + w2));
	  I(i,j,k).dezdx    = fourth*((w3 - w0) + (w1 - w2));
	  I(i,j,k).dezdy    = fourth*((w3 - w0) - (w1 - w2));
	  I(i,j,k).d2ezdxdy = fourth*((w3 + w0) - (w1 + w2));
	  
	  // bx interpolation coefficients
	  w0 = F(i  ,j,k).cbx;
	  w1 = F(i+1,j,k).cbx;
	  I(i,j,k).cbx    = half*(w1 + w0);
	  I(i,j,k).dcbxdx = half*(w1 - w0);
	  
	  // by interpolation coefficients
	  w0 = F(i,j  ,k).cby;
	  w1 = F(i,j+1,k).cby;
	  I(i,j,k).cby    = half*(w1 + w0);
	  I(i,j,k).dcbydy = half*(w1 - w0);
	  
	  // bz interpolation coefficients
	  w0 = F(i,j,k  ).cbz;
	  w1 = F(i,j,k+1).cbz;
	  I(i,j,k).cbz    = half*(w1 + w0);
	  I(i,j,k).dcbzdz = half*(w1 - w0);
	}
      }
    }
  }
  
  void load(FieldArray& fa)
  {
    TIC do_load(fa); TOC(load_interpolator, 1);
  }
};

#endif

