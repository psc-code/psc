
#ifndef PSC_INTERPOLATOR_OPS_H
#define PSC_INTERPOLATOR_OPS_H

template<class IA, class FA>
struct PscInterpolatorOps {
  typedef IA Interpolator;
  typedef FA FieldArray;

  void load_interpolator(Interpolator& ia, /*const*/ FieldArray& fa)
  {
    Field3D<VpicFieldArray> F(fa);
    Field3D<VpicInterpolator> I(ia);

    enum {
      CBX = FieldArray::CBX,
      CBY = FieldArray::CBY,
      CBZ = FieldArray::CBZ,
      EX  = FieldArray::EX,
      EY  = FieldArray::EY,
      EZ  = FieldArray::EZ,
    };
    
    const int nx = ia.g->nx;
    const int ny = ia.g->ny;
    const int nz = ia.g->nz;
    
    const float fourth = 0.25;
    const float half   = 0.5;
    
    float w0, w1, w2, w3;
    
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  
	  // ex interpolation
	  //w0 = f(i,j  ,k  ).ex;
	  w0 = F(EX, i,j  ,k  );
	  w1 = F(EX, i,j+1,k  );
	  w2 = F(EX, i,j  ,k+1);
	  w3 = F(EX, i,j+1,k+1);
	  I(ia, i,j,k).ex       = fourth*((w3 + w0) + (w1 + w2));
	  I(ia, i,j,k).dexdy    = fourth*((w3 - w0) + (w1 - w2));
	  I(ia, i,j,k).dexdz    = fourth*((w3 - w0) - (w1 - w2));
	  I(ia, i,j,k).d2exdydz = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ey interpolation coefficients
	  w0 = F(EY, i  ,j,k  );
	  w1 = F(EY, i  ,j,k+1);
	  w2 = F(EY, i+1,j,k  );
	  w3 = F(EY, i+1,j,k+1);
	  I(ia, i,j,k).ey       = fourth*((w3 + w0) + (w1 + w2));
	  I(ia, i,j,k).deydz    = fourth*((w3 - w0) + (w1 - w2));
	  I(ia, i,j,k).deydx    = fourth*((w3 - w0) - (w1 - w2));
	  I(ia, i,j,k).d2eydzdx = fourth*((w3 + w0) - (w1 + w2));
	  
	  // ez interpolation coefficients
	  w0 = F(EZ, i  ,j  ,k);
	  w1 = F(EZ, i+1,j  ,k);
	  w2 = F(EZ, i  ,j+1,k);
	  w3 = F(EZ, i+1,j+1,k);
	  I(ia, i,j,k).ez       = fourth*((w3 + w0) + (w1 + w2));
	  I(ia, i,j,k).dezdx    = fourth*((w3 - w0) + (w1 - w2));
	  I(ia, i,j,k).dezdy    = fourth*((w3 - w0) - (w1 - w2));
	  I(ia, i,j,k).d2ezdxdy = fourth*((w3 + w0) - (w1 + w2));
	  
	  // bx interpolation coefficients
	  w0 = F(CBX, i  ,j,k);
	  w1 = F(CBX, i+1,j,k);
	  I(ia, i,j,k).cbx    = half*(w1 + w0);
	  I(ia, i,j,k).dcbxdx = half*(w1 - w0);
	  
	  // by interpolation coefficients
	  w0 = F(CBY, i,j  ,k);
	  w1 = F(CBY, i,j+1,k);
	  I(ia, i,j,k).cby    = half*(w1 + w0);
	  I(ia, i,j,k).dcbydy = half*(w1 - w0);
	  
	  // bz interpolation coefficients
	  w0 = F(CBZ, i,j,k  );
	  w1 = F(CBZ, i,j,k+1);
	  I(ia, i,j,k).cbz    = half*(w1 + w0);
	  I(ia, i,j,k).dcbzdz = half*(w1 - w0);
	}
      }
    }
  }
  
  void load_interpolator_array(Interpolator *interpolator,
			       FieldArray *vmflds)
  {
    TIC load_interpolator(*interpolator, *vmflds); TOC(load_interpolator, 1);
  }
};

#endif

