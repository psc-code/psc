
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

  // ----------------------------------------------------------------------
  // synchronize_jf
  
  void synchronize_jf(FieldArray& fa) {
    field_t * f;
    int size, face, x, y, z;
    float *p, lw, rw;

    Field3D<FieldArray> F(fa);
    grid_t* g = fa.g;

    local_adjust_jf(fa.f, g );

    const int nx = fa.g->nx, ny = fa.g->ny, nz = fa.g->nz;

// Generic looping
#define XYZ_LOOP(xl,xh,yl,yh,zl,zh) \
  for( z=zl; z<=zh; z++ )	    \
    for( y=yl; y<=yh; y++ )	    \
      for( x=xl; x<=xh; x++ )
	      
// yz_EDGE_LOOP => Loop over all non-ghost y-oriented edges at plane x
#define yz_EDGE_LOOP(x) XYZ_LOOP(x,x,1,ny,1,nz+1)
#define zx_EDGE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz)
#define xy_EDGE_LOOP(z) XYZ_LOOP(1,nx,1,ny+1,z,z)

// zy_EDGE_LOOP => Loop over all non-ghost z-oriented edges at plane x
#define zy_EDGE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz)
#define xz_EDGE_LOOP(y) XYZ_LOOP(1,nx,y,y,1,nz+1)
#define yx_EDGE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny,z,z)

# define BEGIN_RECV(i,j,k,X,Y,Z)                                        \
    begin_recv_port(i,j,k, ( n##Y*(n##Z+1) +				\
			     n##Z*(n##Y+1) + 1 )*sizeof(float), g )

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {			\
      size = ( n##Y*(n##Z+1) +						\
	       n##Z*(n##Y+1) + 1 )*sizeof(float);			\
	       p = (float *)size_send_port( i, j, k, size, g );		\
	       if( p ) {						\
		 (*(p++)) = g->d##X;					\
		 face = (i+j+k)<0 ? 1 : n##X+1;				\
		 Y##Z##_EDGE_LOOP(face) (*(p++)) = F(x,y,z).jf##Y;	\
		 Z##Y##_EDGE_LOOP(face) (*(p++)) = F(x,y,z).jf##Z;	\
		 begin_send_port( i, j, k, size, g );			\
	       }							\
    } END_PRIMITIVE

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {                \
      p = (float *)end_recv_port(i,j,k,g);			\
      if( p ) {							\
	rw = (*(p++));                 /* Remote g->d##X */	\
	lw = rw + g->d##X;					\
	rw /= lw;						\
	lw = g->d##X/lw;					\
	lw += lw;						\
	rw += rw;						\
	face = (i+j+k)<0 ? n##X+1 : 1; /* Twice weighted sum */	\
	Y##Z##_EDGE_LOOP(face) {				\
	  f = &F(x,y,z);					\
	  f->jf##Y = lw*f->jf##Y + rw*(*(p++));			\
	}							\
	Z##Y##_EDGE_LOOP(face) {				\
	  f = &F(x,y,z);					\
	  f->jf##Z = lw*f->jf##Z + rw*(*(p++));			\
	}							\
      }								\
    } END_PRIMITIVE

# define END_SEND(i,j,k,X,Y,Z) end_send_port( i, j, k, g )

    // Exchange x-faces
    BEGIN_SEND(-1, 0, 0,x,y,z);
    BEGIN_SEND( 1, 0, 0,x,y,z);
    BEGIN_RECV(-1, 0, 0,x,y,z);
    BEGIN_RECV( 1, 0, 0,x,y,z);
    END_RECV(-1, 0, 0,x,y,z);
    END_RECV( 1, 0, 0,x,y,z);
    END_SEND(-1, 0, 0,x,y,z);
    END_SEND( 1, 0, 0,x,y,z);

    // Exchange y-faces
    BEGIN_SEND( 0,-1, 0,y,z,x);
    BEGIN_SEND( 0, 1, 0,y,z,x);
    BEGIN_RECV( 0,-1, 0,y,z,x);
    BEGIN_RECV( 0, 1, 0,y,z,x);
    END_RECV( 0,-1, 0,y,z,x);
    END_RECV( 0, 1, 0,y,z,x);
    END_SEND( 0,-1, 0,y,z,x);
    END_SEND( 0, 1, 0,y,z,x);

    // Exchange z-faces
    BEGIN_SEND( 0, 0,-1,z,x,y);
    BEGIN_SEND( 0, 0, 1,z,x,y);
    BEGIN_RECV( 0, 0,-1,z,x,y);
    BEGIN_RECV( 0, 0, 1,z,x,y);
    END_RECV( 0, 0,-1,z,x,y);
    END_RECV( 0, 0, 1,z,x,y);
    END_SEND( 0, 0,-1,z,x,y);
    END_SEND( 0, 0, 1,z,x,y);

# undef BEGIN_RECV
# undef BEGIN_SEND
# undef END_RECV
# undef END_SEND
  }

  // ----------------------------------------------------------------------
  // compute_div_b_err
  
  static void compute_div_b_err(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    const float px = (nx>1) ? g->rdx : 0;  // FIXME, should be based on global dims
    const float py = (ny>1) ? g->rdy : 0;
    const float pz = (nz>1) ? g->rdz : 0;

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  F(i,j,k).div_b_err = (px * (F(i+1,j,k).cbx - F(i,j,k).cbx) +
				py * (F(i,j+1,k).cby - F(i,j,k).cby) +
				pz * (F(i,j,k+1).cbz - F(i,j,k).cbz));
	}
      }
    }
  }

  // ----------------------------------------------------------------------
  // compute_div_e_err
  
  void compute_div_e_err(FieldArray& fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);

    // Begin setting normal e ghosts
    begin_remote_ghost_norm_e(fa.f, fa.g);
    local_ghost_norm_e(fa.f, fa.g);

    // Overlap local computation
    Field3D<FieldArray> F(fa);
    const material_coefficient_t* m = params->mc;
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    const float nc = m->nonconductive;
    const float px = ((nx>1) ? g->rdx : 0) * m->epsx;
    const float py = ((ny>1) ? g->rdy : 0) * m->epsy;
    const float pz = ((nz>1) ? g->rdz : 0) * m->epsz;
    const float cj = 1./g->eps0;

#define UPDATE_DERR_E(i,j,k) F(i,j,k).div_e_err = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) + \
						      py * (F(i,j,k).ey - F(i,j-1,k).ey) + \
						      pz * (F(i,j,k).ez - F(i,j,k-1).ez) - \
						      cj * (F(i,j,k).rhof + F(i,j,k).rhob) )
    
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  UPDATE_DERR_E(i,j,k);
	}
      }
    }

    // Finish setting normal e ghosts
    end_remote_ghost_norm_e(fa.f, fa.g);

    // z faces, x edges, y edges and all corners
    for(int j = 1; j <= ny+1; j++) {
      for(int i = 1; i <= nx+1; i++) {
	UPDATE_DERR_E(i,j,1   );
      }
    }
    for(int j = 1; j <= ny+1; j++) {
      for(int i = 1; i <= nx+1; i++) {
	UPDATE_DERR_E(i,j,nz+1);
      }
    }

    // y faces, z edges
    for(int k = 2; k <= nz; k++) {
      for(int i = 1; i <= nx+1; i++) {
	UPDATE_DERR_E(i,1   ,k);
      }
    }
    for(int k = 2; k <= nz; k++) {
      for(int i = 1; i <= nx+1; i++) {
	UPDATE_DERR_E(i,ny+1,k);
      }
    }

    // x faces
    for(int k = 2; k <= nz; k++) {
      for(int j = 2; j <= ny; j++) {
	UPDATE_DERR_E(1   ,j,k);
	UPDATE_DERR_E(nx+1,j,k);
      }
    }
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_b_err
  //
  // OPT: doing that at the same time as div_b should be faster

  double compute_rms_div_b_err(FieldArray &fa)
  {
    return fa.kernel->compute_rms_div_b_err(&fa);
  }
  
};


#endif
