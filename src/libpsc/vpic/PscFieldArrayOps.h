
#ifndef PSC_FIELD_ARRAY_OPS_H
#define PSC_FIELD_ARRAY_OPS_H

// ======================================================================
// PscFieldArrayOps

template<class FA>
struct PscFieldArrayOps {
  typedef FA FieldArray;

  // ----------------------------------------------------------------------
  // advance_e
  
  void advance_e(FieldArray& fa, double frac)
  {
    fa.kernel->advance_e(&fa, frac);
  }

  // ----------------------------------------------------------------------
  // advance_b
  
#define CBX FieldArray::CBX
#define CBY FieldArray::CBY
#define CBZ FieldArray::CBZ
#define EX FieldArray::EX
#define EY FieldArray::EY
#define EZ FieldArray::EZ

#define UPDATE_CBX() F(CBX, i,j,k) -= (py*(F(EZ, i,j+1,k) - F(EZ ,i,j,k)) - pz*(F(EY, i,j,k+1) - F(EY, i,j,k)))
#define UPDATE_CBY() F(CBY, i,j,k) -= (pz*(F(EX, i,j,k+1) - F(EX ,i,j,k)) - px*(F(EZ, i+1,j,k) - F(EZ, i,j,k)))
#define UPDATE_CBZ() F(CBZ, i,j,k) -= (px*(F(EY, i+1,j,k) - F(EY, i,j,k)) - py*(F(EX, i,j+1,k) - F(EX, i,j,k)))

  void advance_b(FieldArray& fa, double frac)
  {
    Field3D<FieldArray> F(fa);
    const grid_t *g  = fa.g;
    int nx = g->nx, ny = g->ny, nz = g->nz;

    // FIXME, invariant should be based on global dims
    const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;
    const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;
    const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0;

    // bulk
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBX(); UPDATE_CBY(); UPDATE_CBZ();
	}
      }
    }

    // leftover bx
    { int i = nx + 1;
      for (int k = 1; k <= nz; k++) {
	for (int j = 1; j <= ny; j++) {
	  UPDATE_CBX();
	}
      }
    }
    
    // leftover by
    { int j = ny + 1;
      for (int k = 1; k <= nz; k++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBY();
	}
      }
    }
    
    // leftover bz
    { int k = nz + 1;
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  UPDATE_CBZ();
	}
      }
    }
    
    local_adjust_norm_b(fa.f, fa.g);
  }

#undef CBX
#undef CBY
#undef CBZ
#undef EX
#undef EY
#undef EZ

  // ----------------------------------------------------------------------
  // clear_jf

  void clear_jf(FieldArray& fa)
  {
    const int nv = fa.g->nv;
    for (int v = 0; v < nv; v++) {
      fa[v].jfx = 0;
      fa[v].jfy = 0;
      fa[v].jfz = 0;
    }
  }

  // ----------------------------------------------------------------------
  // clear_rhof

  void clear_rhof(FieldArray& fa)
  {
    const int nv = fa.g->nv;
    for (int v = 0; v < nv; v++) {
      fa[v].rhof = 0;
    }
  }
};


#endif
