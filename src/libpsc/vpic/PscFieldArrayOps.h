
#ifndef PSC_FIELD_ARRAY_OPS_H
#define PSC_FIELD_ARRAY_OPS_H

// ======================================================================
// PscFieldArrayOps

template<class FA>
struct PscFieldArrayOps {
  typedef FA FieldArray;

  // ----------------------------------------------------------------------
  // foreach

  template<class F>
  static void foreach(F f, int ib, int ie, int jb, int je, int kb, int ke)
  {
    for (int k = kb; k <= ke; k++) {
      for (int j = jb; j <= je; j++) {
	for (int i = ib; i <= ie; i++) {
	  f(i,j,k);
	}
      }
    }
  };

  template<class F>
  static void foreach_nc_interior(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    foreach(f, 2, nx, 2, ny, 2, nz);
  }

  template<class F>
  static void foreach_nc_boundary(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    // z faces, x edges, y edges and all corners
    foreach(f, 1   , nx+1, 1   , ny+1, 1   , 1   );
    foreach(f, 1   , nx+1, 1   , ny+1, nz+1, nz+1);

    // y faces, z edges
    foreach(f, 1   , nx+1, 1   , 1   ,    2, nz  );
    foreach(f, 1   , nx+1, ny+1, ny+1,    2, nz  );

    // x faces
    foreach(f, 1   , 1   , 2   , ny  ,    2, nz  );
    foreach(f, nx+1, nx+1, 2   , ny  ,    2, nz  );
  }

  template<class F>
  static void foreach_ec_interior(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  f.x(i,j,k);
	  f.y(i,j,k);
	  f.z(i,j,k);
	}
      }
    }

    // Do leftover interior ex
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.x(1,j,k);
      }
    }

    // Do leftover interior ey
    for (int k = 2; k <= nz; k++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,1,k);
      }
    }

    // Do leftover interior ez
    for (int j = 2; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.z(i,j,1);
      }
    }
  }

  template<class F>
  static void foreach_ec_boundary(F f, const grid_t *g)
  {
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    for (int j = 1; j <= ny+1; j++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,j,1);
      }
    }
    for (int j = 1; j <= ny+1; j++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,j,nz+1);
      }
    }
    for (int k = 2; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,1,k);
      }
    }
    for (int k = 2; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
	f.x(i,ny+1,k);
      }
    }

    // Do exterior ey
    for (int k = 1; k <= nz+1; k++) {
      for (int j = 1; j <= ny; j++) {
	f.y(1,j,k);
      }
    }
    for (int k = 1; k <= nz+1; k++) {
      for (int j = 1; j <= ny; j++) {
	f.y(nx+1,j,k);
      }
    }
    for (int j = 1; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,j,1);
      }
    }
    for (int j = 1; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.y(i,j,nz+1);
      }
    }

    // Do exterior ez
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx+1; i++) {
	f.z(i,1,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx+1; i++) {
	f.z(i,ny+1,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.z(1,j,k);
      }
    }
    for (int k = 1; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	f.z(nx+1,j,k);
      }
    }
  }
  
  // ----------------------------------------------------------------------
  // advance_e
  
  void vacuum_advance_e(FieldArray& fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);

    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    const material_coefficient_t* m = params->mc;

    struct AdvanceE {
      AdvanceE(FieldArray& fa, const grid_t *g, const material_coefficient_t *m,
	       const double damp_)
	: F(fa),
	  decayx(m->decayx),
	  decayy(m->decayy),
	  decayz(m->decayz),
	  drivex(m->drivex),
	  drivey(m->drivey),
	  drivez(m->drivez),
	  px_muz(g->nx > 1 ? (1+damp_)*g->cvac*g->dt*g->rdx*m->rmuz : 0),
	  px_muy(g->nx > 1 ? (1+damp_)*g->cvac*g->dt*g->rdx*m->rmuy : 0),
	  py_mux(g->ny > 1 ? (1+damp_)*g->cvac*g->dt*g->rdy*m->rmux : 0),
	  py_muz(g->ny > 1 ? (1+damp_)*g->cvac*g->dt*g->rdy*m->rmuz : 0),
	  pz_muy(g->nz > 1 ? (1+damp_)*g->cvac*g->dt*g->rdz*m->rmuy : 0),
	  pz_mux(g->nz > 1 ? (1+damp_)*g->cvac*g->dt*g->rdz*m->rmux : 0),
	  damp(damp_),
	  cj(g->dt / g->eps0)
      {
      }

      void x(int i, int j, int k)
      {
	F(i,j,k).tcax = ((py_muz*(F(i,j,k).cbz - F(i,j-1,k).cbz) -
			  pz_muy*(F(i,j,k).cby - F(i,j,k-1).cby)) -
			 damp*F(i,j,k).tcax);
	F(i,j,k).ex   = decayx*F(i,j,k).ex + drivex*(F(i,j,k).tcax - cj*F(i,j,k).jfx);
      }
      
      void y(int i, int j, int k)
      {
	F(i,j,k).tcay = ((pz_mux*(F(i,j,k).cbx - F(i,j,k-1).cbx) -
			  px_muz*(F(i,j,k).cbz - F(i-1,j,k).cbz)) -
			 damp*F(i,j,k).tcay);
	F(i,j,k).ey   = decayy*F(i,j,k).ey + drivey*(F(i,j,k).tcay - cj*F(i,j,k).jfy);
      }
      
      void z(int i, int j, int k)
      {
	F(i,j,k).tcaz = ((px_muy*(F(i,j,k).cby - F(i-1,j,k).cby) -
			  py_mux*(F(i,j,k).cbx - F(i,j-1,k).cbx)) -
			 damp*F(i,j,k).tcaz);
	F(i,j,k).ez   = decayz*F(i,j,k).ez + drivez*(F(i,j,k).tcaz - cj*F(i,j,k).jfz);
      }

      Field3D<FieldArray> F;
      const float decayx, decayy, decayz, drivex, drivey, drivez;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
      const float damp, cj;
    };

    AdvanceE advanceE(fa, fa.g, m, params->damp);

    begin_remote_ghost_tang_b(fa.f, fa.g);

    local_ghost_tang_b(fa.f, fa.g);
    foreach_ec_interior(advanceE, fa.g);

    end_remote_ghost_tang_b(fa.f, fa.g);

    foreach_ec_boundary(advanceE, fa.g);
    local_adjust_tang_e(fa.f, fa.g);
  }

  void advance_e(FieldArray& fa, double frac)
  {
    assert(frac == 1.);
    vacuum_advance_e(fa);
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

// x_NODE_LOOP => Loop over all non-ghost nodes at plane x
#define x_NODE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz+1)
#define y_NODE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz+1)
#define z_NODE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny+1,z,z)

// x_FACE_LOOP => Loop over all x-faces at plane x
#define x_FACE_LOOP(x) XYZ_LOOP(x,x,1,ny,1,nz)
#define y_FACE_LOOP(y) XYZ_LOOP(1,nx,y,y,1,nz)
#define z_FACE_LOOP(z) XYZ_LOOP(1,nx,1,ny,z,z)

  // ----------------------------------------------------------------------
  // synchronize_tang_e_norm_b

  double synchronize_tang_e_norm_b(FieldArray& fa)
  {
    float * p;
    double w1, w2, err = 0, gerr;
    int size, face, x, y, z;

    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    
    local_adjust_tang_e(fa.f, fa.g);
    local_adjust_norm_b(fa.f, fa.g);

# define BEGIN_RECV(i,j,k,X,Y,Z)					\
    begin_recv_port(i,j,k, ( 2*n##Y*(n##Z+1) + 2*n##Z*(n##Y+1) +	\
			     n##Y*n##Z )*sizeof(float), g )

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {              \
      size = ( 2*n##Y*(n##Z+1) + 2*n##Z*(n##Y+1) +		\
	       n##Y*n##Z )*sizeof(float);			\
      p = (float *)size_send_port( i, j, k, size, g );		\
      if( p ) {							\
	face = (i+j+k)<0 ? 1 : n##X+1;				\
	X##_FACE_LOOP(face) (*(p++)) = F(x,y,z).cb##X;	\
	Y##Z##_EDGE_LOOP(face) {				\
	  field_t* f = &F(x,y,z);				\
	  (*(p++)) = f->e##Y;					\
	  (*(p++)) = f->tca##Y;					\
	}							\
	Z##Y##_EDGE_LOOP(face) {				\
	  field_t* f = &F(x,y,z);				\
	  (*(p++)) = f->e##Z;					\
	  (*(p++)) = f->tca##Z;					\
	}							\
	begin_send_port( i, j, k, size, g );			\
      }								\
    } END_PRIMITIVE
    
# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {		\
      p = (float *)end_recv_port(i,j,k,g);			\
      if( p ) {							\
	face = (i+j+k)<0 ? n##X+1 : 1; /* Average */		\
	X##_FACE_LOOP(face) {					\
	  field_t* f = &F(x,y,z);				\
	  w1 = (*(p++));					\
	  w2 = f->cb##X;					\
	  f->cb##X = 0.5*( w1+w2 );				\
	  err += (w1-w2)*(w1-w2);				\
	}							\
	Y##Z##_EDGE_LOOP(face) {				\
	  field_t* f = &F(x,y,z);				\
	  w1 = (*(p++));					\
	  w2 = f->e##Y;						\
	  f->e##Y = 0.5*( w1+w2 );				\
	  err += (w1-w2)*(w1-w2);				\
	  w1 = (*(p++));					\
	  w2 = f->tca##Y;					\
	  f->tca##Y = 0.5*( w1+w2 );				\
	}							\
	Z##Y##_EDGE_LOOP(face) {				\
	  field_t* f = &F(x,y,z);				\
	  w1 = (*(p++));					\
	  w2 = f->e##Z;						\
	  f->e##Z = 0.5*( w1+w2 );				\
	  err += (w1-w2)*(w1-w2);				\
	  w1 = (*(p++));					\
	  w2 = f->tca##Z;					\
	  f->tca##Z = 0.5*( w1+w2 );				\
	}							\
      }								\
    } END_PRIMITIVE

# define END_SEND(i,j,k,X,Y,Z) end_send_port( i, j, k, g )

    // Exchange x-faces
    BEGIN_RECV(-1, 0, 0,x,y,z);
    BEGIN_RECV( 1, 0, 0,x,y,z);
    BEGIN_SEND(-1, 0, 0,x,y,z);
    BEGIN_SEND( 1, 0, 0,x,y,z);
    END_SEND(-1, 0, 0,x,y,z);
    END_SEND( 1, 0, 0,x,y,z);
    END_RECV(-1, 0, 0,x,y,z);
    END_RECV( 1, 0, 0,x,y,z);

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

    mp_allsum_d( &err, &gerr, 1 );
    return gerr;
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
  // synchronize_rho
  
  // Note: synchronize_rho assumes that rhof has _not_ been adjusted at
  // the local domain boundary to account for partial cells but that
  // rhob _has_.  Specifically it is very expensive to accumulate rhof
  // and doing the adjustment for each particle is adds even more
  // expense.  Worse, if we locally corrected it after each species,
  // we cannot accumulate the next species in the same unless we use
  // (running sum of locally corrected results and thw current species
  // rhof being accumulated).  Further, rhof is always accumulated from
  // scratch so we don't have to worry about whether or not the previous
  // rhof values were in a locally corrected form.  Thus, after all
  // particles have accumulated to rhof, we correct it for partial cells
  // and remote cells for use with divergence cleaning and this is
  // the function that does the correction.
  //
  // rhob is another story though.  rhob is continuously incrementally
  // accumulated over time typically through infrequent surface area
  // scaling processes.  Like rho_f, after synchronize_rhob, rhob _must_
  // be corrected for partial and remote celle for the benefit of
  // divergence cleaning. And like rho_f, since we don't want to have
  // to keep around two versions of rhob (rhob contributions since last
  // sync and rhob as of last sync), we have no choice but to do the
  // charge accumulation per particle to rhob in a locally corrected
  // form.

  void synchronize_rho(FieldArray& fa) {
    int size, face, x, y, z;
    float *p, hlw, hrw, lw, rw;

    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    local_adjust_rhof(fa.f, fa.g);
    local_adjust_rhob(fa.f, fa.g);

# define BEGIN_RECV(i,j,k,X,Y,Z)					\
    begin_recv_port(i,j,k, ( 1 + 2*(n##Y+1)*(n##Z+1) )*sizeof(float), g )

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {      \
      size = ( 1 + 2*(n##Y+1)*(n##Z+1) )*sizeof(float);	\
      p = (float *)size_send_port( i, j, k, size, g );	\
      if( p ) {						\
	(*(p++)) = g->d##X;				\
	face = (i+j+k)<0 ? 1 : n##X+1;			\
	X##_NODE_LOOP(face) {				\
	  field_t *f = &F(x,y,z);			\
	  (*(p++)) = f->rhof;				\
	  (*(p++)) = f->rhob;				\
	}						\
	begin_send_port( i, j, k, size, g );		\
      }							\
    } END_PRIMITIVE

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {                \
      p = (float *)end_recv_port(i,j,k,g);			\
      if( p ) {							\
	hrw  = (*(p++));               /* Remote g->d##X */	\
	hlw  = hrw + g->d##X;					\
	hrw /= hlw;						\
	hlw  = g->d##X/hlw;					\
	lw   = hlw + hlw;					\
	rw   = hrw + hrw;					\
	face = (i+j+k)<0 ? n##X+1 : 1;				\
	X##_NODE_LOOP(face) {					\
	  field_t *f = &F(x,y,z);				\
	  f->rhof =  lw*f->rhof  + rw*(*(p++));			\
	  f->rhob = hlw*f->rhob + hrw*(*(p++));			\
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
  // compute_rhob
  void compute_rhob(FieldArray& fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);
    const material_coefficient_t* m = params->mc;

    struct CalcRhoB {
      CalcRhoB(FieldArray& fa, const material_coefficient_t *m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.g->nx > 1 ? fa.g->eps0 * m->epsx * fa.g->rdx : 0),
	  py(fa.g->ny > 1 ? fa.g->eps0 * m->epsy * fa.g->rdy : 0),
	  pz(fa.g->nz > 1 ? fa.g->eps0 * m->epsz * fa.g->rdz : 0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).rhob = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
			    py * (F(i,j,k).ey - F(i,j-1,k).ey) +
			    pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
			    F(i,j,k).rhof);
      }

      Field3D<FieldArray> F;
      const float nc, px, py, pz;
    };

    CalcRhoB updater(fa, m);

    // Begin setting normal e ghosts
    begin_remote_ghost_norm_e(fa.f, fa.g);

    // Overlap local computation
    local_ghost_norm_e(fa.f, fa.g);
    foreach_nc_interior(updater, fa.g);
    
    // Finish setting normal e ghosts
    end_remote_ghost_norm_e(fa.f, fa.g);

    // Now do points on boundary
    foreach_nc_boundary(updater, fa.g);

    local_adjust_rhob(fa.f, fa.g);
  }

  // ----------------------------------------------------------------------
  // compute_div_e_err
  
  void compute_div_e_err(FieldArray& fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);
    const material_coefficient_t* m = params->mc;

    struct CalcDivE {
      CalcDivE(FieldArray& fa, const material_coefficient_t *m)
	: F(fa),
	  nc(m->nonconductive),
	  px(fa.g->nx > 1 ? fa.g->eps0 * fa.g->rdx : 0),
	  py(fa.g->ny > 1 ? fa.g->eps0 * fa.g->rdy : 0),
	  pz(fa.g->nz > 1 ? fa.g->eps0 * fa.g->rdz : 0),
	  cj(1. / fa.g->eps0)
      {
      }

      void operator()(int i, int j, int k)
      {
	F(i,j,k).div_e_err = nc*(px * (F(i,j,k).ex - F(i-1,j,k).ex) +
				 py * (F(i,j,k).ey - F(i,j-1,k).ey) +
				 pz * (F(i,j,k).ez - F(i,j,k-1).ez) -
				 cj * (F(i,j,k).rhof + F(i,j,k).rhob));
      }

      Field3D<FieldArray> F;
      const float nc, px, py, pz, cj;
    };

    CalcDivE updater(fa, m);
    
    // Begin setting normal e ghosts
    begin_remote_ghost_norm_e(fa.f, fa.g);

    // Overlap local computation
    local_ghost_norm_e(fa.f, fa.g);
    foreach_nc_interior(updater, fa.g);

    // Finish setting normal e ghosts
    end_remote_ghost_norm_e(fa.f, fa.g);

    // Now do points on boundary
    foreach_nc_boundary(updater, fa.g);

    local_adjust_div_e(fa.f, fa.g);
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_b_err
  //
  // OPT: doing that at the same time as div_b should be faster

  double compute_rms_div_b_err(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    const int nx = fa.g->nx, ny = fa.g->ny, nz = fa.g->nz;

    double err = 0;
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  err += sqr(F(i,j,k).div_b_err);
	}
      }
    }
    
    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * fa.g->dV;
    local[1] = (nx * ny * nz) * fa.g->dV;
    mp_allsum_d( local, _global, 2 );
    return fa.g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // compute_rms_div_e_err
  
  double compute_rms_div_e_err(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    const int nx = fa.g->nx, ny = fa.g->ny, nz = fa.g->nz;

    // Interior points
    // FIXME, it's inconsistent to calc the square in single prec here, but
    // double prec later
    double err = 0;
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  err += sqr(F(i,j,k).div_e_err);
	}
      }
    }

    int x, y, z;
    // Exterior faces

    for( y=2; y<=ny; y++ ) {
      for( z=2; z<=nz; z++ ) {
	err += 0.5*sqr((double) (F(1   ,y,z).div_e_err));
	err += 0.5*sqr((double) (F(nx+1,y,z).div_e_err));
      }
    }

    for( z=2; z<=nz; z++ ) {
      for( x=2; x<=nx; x++ ) {
	err += 0.5*sqr((double) (F(x,1   ,z).div_e_err));
	err += 0.5*sqr((double) (F(x,ny+1,z).div_e_err));
      }
    }
    
    for( x=2; x<=nx; x++ ) {
      for( y=2; y<=ny; y++ ) {
	err += 0.5*sqr((double) (F(x,y,1   ).div_e_err));
	err += 0.5*sqr((double) (F(x,y,nz+1).div_e_err));
      }
    }

    // Exterior edges

    for( x=2; x<=nx; x++ ) {
      err += 0.25*sqr((double) (F(x,1   ,1   ).div_e_err));
      err += 0.25*sqr((double) (F(x,ny+1,1   ).div_e_err));
      err += 0.25*sqr((double) (F(x,1   ,nz+1).div_e_err));
      err += 0.25*sqr((double) (F(x,ny+1,nz+1).div_e_err));
    }

    for( y=2; y<=ny; y++ ) {
      err += 0.25*sqr((double) (F(1   ,y,1   ).div_e_err));
      err += 0.25*sqr((double) (F(1   ,y,nz+1).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,y,1   ).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,y,nz+1).div_e_err));
    }

    for( z=2; z<=nz; z++ ) {
      err += 0.25*sqr((double) (F(1   ,1   ,z).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,1   ,z).div_e_err));
      err += 0.25*sqr((double) (F(1   ,ny+1,z).div_e_err));
      err += 0.25*sqr((double) (F(nx+1,ny+1,z).div_e_err));
    }

    // Exterior corners

    err += 0.125*sqr((double) (F(1   ,1   ,   1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,1   ,   1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,ny+1,   1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,ny+1,   1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,1   ,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,1   ,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(1   ,ny+1,nz+1).div_e_err));
    err += 0.125*sqr((double) (F(nx+1,ny+1,nz+1).div_e_err));
    
    double local[2], _global[2]; // FIXME, name clash with global macro
    local[0] = err * fa.g->dV;
    local[1] = (nx * ny * nz) * fa.g->dV;
    mp_allsum_d( local, _global, 2 );
    return fa.g->eps0 * sqrt(_global[0]/_global[1]);
  }

  // ----------------------------------------------------------------------
  // clean_div_b
  
#define MARDER_CBX(i,j,k) F(i,j,k).cbx += px * (F(i,j,k).div_b_err - F(i-1,j,k).div_b_err)
#define MARDER_CBY(i,j,k) F(i,j,k).cby += py * (F(i,j,k).div_b_err - F(i,j-1,k).div_b_err)
#define MARDER_CBZ(i,j,k) F(i,j,k).cbz += pz * (F(i,j,k).div_b_err - F(i,j,k-1).div_b_err)

  void clean_div_b(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;

    const int nx = g->nx, ny = g->ny, nz = g->nz;
    float px = (nx>1) ? g->rdx : 0;
    float py = (ny>1) ? g->rdy : 0;
    float pz = (nz>1) ? g->rdz : 0;
    float alphadt = 0.3888889/(px*px + py*py + pz*pz);
    px *= alphadt;
    py *= alphadt;
    pz *= alphadt;

    // Have pipelines do Marder pass in interior.  The host handles
    // stragglers.

    // Begin setting derr ghosts
    begin_remote_ghost_div_b(fa.f, g);
    local_ghost_div_b(fa.f, g);

    // Interior
    for (int k = 2; k <= nz; k++) {
      for (int j = 2; j <= ny; j++) {
	for (int i = 2; i <= nx; i++) {
	  MARDER_CBX(i,j,k);
	  MARDER_CBY(i,j,k);
	  MARDER_CBZ(i,j,k);
	}
      }
    }
  
    int x, y, z;

    // Do left over interior bx
    for(y=1; y<=ny; y++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBX(x,y,1);
      }
    }
    for(z=2; z<=nz; z++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBX(x,1,z);
      }
    }

    // Left over interior by
    for(z=1; z<=nz; z++) {
      for(y=2; y<=ny; y++) {
	MARDER_CBY(1,y,z);
      }
    }
    for(y=2; y<=ny; y++) {
      for(x=2; x<=nx; x++) {
	MARDER_CBY(x,y,1);
      }
    }

    // Left over interior bz
    for(z=2; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,1,z);
      }
    }
    for(z=2; z<=nz; z++) {
      for(y=2; y<=ny; y++) {
	MARDER_CBZ(1,y,z);
      }
    }

    // Finish setting derr ghosts
  
    end_remote_ghost_div_b(fa.f, g);

    // Do Marder pass in exterior

    // Exterior bx
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_CBX(1,y,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_CBX(nx+1,y,z);
      }
    }

    // Exterior by
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBY(x,1,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBY(x,ny+1,z);
      }
    }

    // Exterior bz
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,y,1);
      }
    }
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_CBZ(x,y,nz+1);
      }
    }

    local_adjust_norm_b(fa.f, g);
  }

  // ----------------------------------------------------------------------
  // clean_div_e
  
#define MARDER_EX(i,j,k) F(i,j,k).ex += px * (F(i+1,j,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EY(i,j,k) F(i,j,k).ey += py * (F(i,j+1,k).div_e_err - F(i,j,k).div_e_err)
#define MARDER_EZ(i,j,k) F(i,j,k).ez += pz * (F(i,j,k+1).div_e_err - F(i,j,k).div_e_err)

  void
  vacuum_clean_div_e(FieldArray &fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);

    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;

    const float _rdx = (nx>1) ? g->rdx : 0;
    const float _rdy = (ny>1) ? g->rdy : 0;
    const float _rdz = (nz>1) ? g->rdz : 0;
    const float alphadt = 0.3888889/(_rdx*_rdx + _rdy*_rdy + _rdz*_rdz);
    const float px = (alphadt*_rdx) * params->mc[0].drivex;
    const float py = (alphadt*_rdy) * params->mc[0].drivey;
    const float pz = (alphadt*_rdz) * params->mc[0].drivez;
                     
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	for (int i = 1; i <= nx; i++) {
	  MARDER_EX(i,j,k);
	  MARDER_EY(i,j,k);
	  MARDER_EZ(i,j,k);
	}
      }
    }
    
    int x, y, z;
  
    // Do left over ex
    for(y=1; y<=ny+1; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_EX(x,y,nz+1);
      }
    }
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx; x++) {
	MARDER_EX(x,ny+1,z);
      }
    }

    // Do left over ey
    for(z=1; z<=nz+1; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_EY(nx+1,y,z);
      }
    }
    for(y=1; y<=ny; y++) {
      for(x=1; x<=nx; x++) {
	MARDER_EY(x,y,nz+1);
      }
    }

    // Do left over ez
    for(z=1; z<=nz; z++) {
      for(x=1; x<=nx+1; x++) {
	MARDER_EZ(x,ny+1,z);
      }
    }
    for(z=1; z<=nz; z++) {
      for(y=1; y<=ny; y++) {
	MARDER_EZ(nx+1,y,z);
      }
    }

    local_adjust_tang_e(fa.f, fa.g);
  }

  void clean_div_e(FieldArray& fa)
  {
    vacuum_clean_div_e(fa);
  }

  // ----------------------------------------------------------------------
  // compute_curl_b
  
  void vacuum_compute_curl_b(FieldArray& fa)
  {
    sfa_params_t* params = static_cast<sfa_params_t*>(fa.params);
    assert(params->n_mc == 1);

    // Update interior fields
    // Note: ex all (1:nx,  1:ny+1,1,nz+1) interior (1:nx,2:ny,2:nz)
    // Note: ey all (1:nx+1,1:ny,  1:nz+1) interior (2:nx,1:ny,2:nz)
    // Note: ez all (1:nx+1,1:ny+1,1:nz ) interior (1:nx,1:ny,2:nz)

    const material_coefficient_t* m = params->mc;

    struct CurlB {
      CurlB(FieldArray& fa, const grid_t *g, const material_coefficient_t *m)
	: F(fa),
	  px_muz(g->nx > 1 ? g->cvac*g->dt*g->rdx*m->rmuz : 0),
	  px_muy(g->nx > 1 ? g->cvac*g->dt*g->rdx*m->rmuy : 0),
	  py_mux(g->ny > 1 ? g->cvac*g->dt*g->rdy*m->rmux : 0),
	  py_muz(g->ny > 1 ? g->cvac*g->dt*g->rdy*m->rmuz : 0),
	  pz_muy(g->nz > 1 ? g->cvac*g->dt*g->rdz*m->rmuy : 0),
	  pz_mux(g->nz > 1 ? g->cvac*g->dt*g->rdz*m->rmux : 0)
      {
      }

      void x(int i, int j, int k)
      {
	F(i,j,k).tcax = (py_muz*(F(i,j,k).cbz - F(i,j-1,k).cbz) -
			 pz_muy*(F(i,j,k).cby - F(i,j,k-1).cby));
      }
      
      void y(int i, int j, int k)
      {
	F(i,j,k).tcay = (pz_mux*(F(i,j,k).cbx - F(i,j,k-1).cbx) -
			 px_muz*(F(i,j,k).cbz - F(i-1,j,k).cbz));
      }
      
      void z(int i, int j, int k)
      {
	F(i,j,k).tcaz = (px_muy*(F(i,j,k).cby - F(i-1,j,k).cby) -
			 py_mux*(F(i,j,k).cbx - F(i,j-1,k).cbx));
      }

      Field3D<FieldArray> F;
      const float px_muz, px_muy, py_mux, py_muz, pz_muy, pz_mux;
    };

    CurlB curlB(fa, fa.g, m);
      
    begin_remote_ghost_tang_b(fa.f, fa.g);

    local_ghost_tang_b(fa.f, fa.g);
    foreach_ec_interior(curlB, fa.g);

    end_remote_ghost_tang_b(fa.f, fa.g);

    foreach_ec_boundary(curlB, fa.g);
    local_adjust_tang_e(fa.f, fa.g); // FIXME, is this right here?
  }

  void compute_curl_b(FieldArray& fa)
  {
    vacuum_compute_curl_b(fa);
  }

};


#endif
