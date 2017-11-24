
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
  
#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

  typedef struct pipeline_args {
    field_t      * ALIGNED(128) f;
    const grid_t *              g;
  } pipeline_args_t;

  static void compute_div_b_err_pipeline( pipeline_args_t * args,
					  int pipeline_rank,
					  int n_pipeline )
  {
    field_t      * ALIGNED(128) f = args->f;
    const grid_t *              g = args->g;
  
    field_t * ALIGNED(16) f0;
    field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
    int x, y, z, n_voxel;

    const int nx = g->nx;
    const int ny = g->ny;
    const int nz = g->nz;

    const float px = (nx>1) ? g->rdx : 0;
    const float py = (ny>1) ? g->rdy : 0;
    const float pz = (nz>1) ? g->rdz : 0;

    // Process the voxels assigned to this pipeline
  
    DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
		       pipeline_rank, n_pipeline,
		       x, y, z, n_voxel );

# define LOAD_STENCIL()				\
    f0 = &f(x,  y,  z  );			\
    fx = &f(x+1,y,  z  );			\
    fy = &f(x,  y+1,z  );			\
    fz = &f(x,  y,  z+1)

    LOAD_STENCIL();

    for( ; n_voxel; n_voxel-- ) {
      f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	py*( fy->cby - f0->cby ) +
	pz*( fz->cbz - f0->cbz );
      f0++; fx++; fy++; fz++;
    
      x++;
      if( x>nx ) {
	x=1, y++;
	if( y>ny ) y=1, z++;
	LOAD_STENCIL();
      }
    }

# undef LOAD_STENCIL

  }

#define HAS_V4_PIPELINE
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  static void compute_div_b_err_pipeline_v4(pipeline_args_t * args,
					    int pipeline_rank,
					    int n_pipeline)
  {
    using namespace v4;

    field_t      * ALIGNED(128) f = args->f;
    const grid_t *              g = args->g;

    field_t * ALIGNED(16) f0;
    field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
    int x, y, z, n_voxel;

    const int nx = g->nx;
    const int ny = g->ny;
    const int nz = g->nz;

    const float px = (nx>1) ? g->rdx : 0;
    const float py = (ny>1) ? g->rdy : 0;
    const float pz = (nz>1) ? g->rdz : 0;

    const v4float vpx(px);
    const v4float vpy(py);
    const v4float vpz(pz);

    v4float f0_cbx, f0_cby, f0_cbz; // Voxel quad magnetic fields
    v4float f0_div_b_err;           // Voxel quad div b errs
    v4float fx_cbx;                 // Voxel quad +x neighbor x magnetic fields
    v4float fy_cby;                 // Voxel quad +y neighbor y magnetic fields
    v4float fz_cbz;                 // Voxel quad +z neighbor z magnetic fields

    field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad
    field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors
    field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +x neighbors
    field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +x neighbors

    // Process the voxels assigned to this pipeline 
  
    DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
		       pipeline_rank, n_pipeline,
		       x, y, z, n_voxel );

    // Process bulk of voxels 4 at a time

# define LOAD_STENCIL()				\
    f0 = &f(x,  y,  z  );			\
    fx = &f(x+1,y,  z  );			\
    fy = &f(x,  y+1,z  );			\
    fz = &f(x,  y,  z+1)

# define NEXT_STENCIL(n)			\
    f0##n = f0++;				\
    fx##n = fx++;				\
    fy##n = fy++;				\
    fz##n = fz++;				\
    x++;					\
    if( x>nx ) {				\
      x=1, y++;					\
      if( y>ny ) y=1, z++;			\
      LOAD_STENCIL();				\
    }

    LOAD_STENCIL();

    for( ; n_voxel>3; n_voxel-=4 ) {
      NEXT_STENCIL(0); NEXT_STENCIL(1); NEXT_STENCIL(2); NEXT_STENCIL(3);

      load_4x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx, f0_cbx, f0_cby, f0_cbz );

      fx_cbx = v4float( fx0->cbx, fx1->cbx, fx2->cbx, fx3->cbx );
      fy_cby = v4float( fy0->cby, fy1->cby, fy2->cby, fy3->cby );
      fz_cbz = v4float( fz0->cbz, fz1->cbz, fz2->cbz, fz3->cbz );

      f0_div_b_err = fma( vpx,fx_cbx-f0_cbx, fma( vpy,fy_cby-f0_cby, vpz*(fz_cbz-f0_cbz) ) );

      store_4x1_tr( f0_div_b_err, &f00->div_b_err, &f01->div_b_err, &f02->div_b_err, &f03->div_b_err );
    }

# undef NEXT_STENCIL
# undef LOAD_STENCIL

  }

#endif

  void
  compute_div_b_err( field_array_t * RESTRICT fa ) {
    pipeline_args_t args[1];

    if( !fa ) ERROR(( "Bad args" ));
  
# if 0 // Original non-pipelined version
    for( z=1; z<=nz; z++ ) {
      for( y=1; y<=ny; y++ ) {
	f0 = &f(1,y,z);
	fx = &f(2,y,z);
	fy = &f(1,y+1,z);
	fz = &f(1,y,z+1);
	for( x=1; x<=nx; x++ ) {
	  f0->div_b_err = px*( fx->cbx - f0->cbx ) +
	    py*( fy->cby - f0->cby ) +
	    pz*( fz->cbz - f0->cbz );
	  f0++; fx++; fy++; fz++;
	}
      }
    }
# endif

    args->f = fa->f;
    args->g = fa->g;
    
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)
    compute_div_b_err_pipeline_v4(args, 0, 1);
#else
    compute_div_b_err_pipeline(args, 0, 1);
#endif
    compute_div_b_err_pipeline(args, 1, 1);
    /* EXEC_PIPELINES( compute_div_b_err, args, 0 ); */
    /* WAIT_PIPELINES(); */
  }

  void compute_div_b_err(FieldArray& fa)
  {
    compute_div_b_err(&fa);
  }
  
};


#endif
