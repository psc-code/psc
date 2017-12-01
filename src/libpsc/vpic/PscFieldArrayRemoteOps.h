
#ifndef PSC_FIELD_ARRAY_REMOTE_OPS_H
#define PSC_FIELD_ARRAY_REMOTE_OPS_H

// Indexing macros
#define field(x,y,z) field[ VOXEL(x,y,z, nx,ny,nz) ]

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

// ======================================================================
// PscFieldArrayRemoteOps

template<class FA>
struct PscFieldArrayRemoteOps {
  typedef FA FieldArray;

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
  static void foreach_edge(const grid_t *g, int Y, int Z, int face, F f)
  {
    if        (Y == 1 && Z == 2) { foreach(f, face, face, 1, g->ny  , 1, g->nz+1);
    } else if (Y == 2 && Z == 1) { foreach(f, face, face, 1, g->ny+1, 1, g->nz  );
    } else if (Y == 2 && Z == 0) { foreach(f, 1, g->nx+1, face, face, 1, g->nz  );
    } else if (Y == 0 && Z == 2) { foreach(f, 1, g->nx  , face, face, 1, g->nz+1);
    } else if (Y == 0 && Z == 1) { foreach(f, 1, g->nx  , 1, g->ny+1, face, face);
    } else if (Y == 1 && Z == 0) { foreach(f, 1, g->nx+1, 1, g->ny  , face, face);
    } else {
      assert(0);
    }
  }

 template<class F>
  static void foreach_node(const grid_t *g, int X, int face, F f)
  {
    if        (X == 0) { foreach(f, face, face, 1, g->ny+1, 1, g->nz+1);
    } else if (X == 1) { foreach(f, 1, g->nx+1, face, face, 1, g->nz+1);
    } else if (X == 2) { foreach(f, 1, g->nx+1, 1, g->ny+1, face, face);
    } else {
      assert(0);
    }
  }

  // ----------------------------------------------------------------------
  // Comm

  struct Comm
  {
    Comm(const grid_t* g) : g_(g)
    {
      nx_[0] = g_->nx;
      nx_[1] = g_->ny;
      nx_[2] = g_->nz;
    }

    // wrap what probably should be in Grid::

    float* get_send_buf(int i, int j, int k, int sz) const
    {
      return static_cast<float*>(::size_send_port(i, j, k, sz * sizeof(float), g_));
    }

    void begin_send_port(int i, int j, int k, int sz) const
    {
      ::begin_send_port(i, j, k, sz * sizeof(float), g_);
    }

    void begin_recv_port(int i, int j, int k, int sz) const
    {
      ::begin_recv_port(i, j, k, sz * sizeof(float), g_);
    }

    void end_send(int i, int j, int k) const
    {
      ::end_send_port(i, j, k, g_);
    }

    // ghost communication
    
    void nc_begin_recv(int i, int j, int k, int X, int Y, int Z) const
    {
      int sz = 1 + (nx_[Y] + 1) * (nx_[Z] + 1);
      begin_recv_port(i, j, k, sz);
    }

    void ec_begin_recv(int i, int j, int k, int X, int Y, int Z) const
    {
      int sz = 1 + nx_[Y] * (nx_[Z]+1) + nx_[Z] * (nx_[Y]+1);
      begin_recv_port(i, j, k, sz);
    }

    template<class F3D>
    void nc_begin_send(int i, int j, int k, int X, int Y, int Z, F3D& F) const
    {
      int sz = 1 + (nx_[Y] + 1) * (nx_[Z] + 1);
      float *p = get_send_buf(i, j, k, sz);
      if (p) {
	int face = (i+j+k) < 0 ? 1 : nx_[X];
	*p++ = (&g_->dx)[X];
	foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).ex)[X]; });
	begin_send_port(i, j, k, sz);
      }
    }
    
    template<class F3D>
    void ec_begin_send(int i, int j, int k, int X, int Y, int Z, F3D& F) const
    {
      int sz = 1 + nx_[Y] * (nx_[Z]+1) + nx_[Z] * (nx_[Y]+1);
      float *p = get_send_buf(i, j, k, sz);
      if (p) {
	int face = (i+j+k) < 0 ? 1 : nx_[X];
	*p++ = (&g_->dx)[X];
	foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Y]; });
	foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Z]; });
	begin_send_port(i, j, k, sz);
      }
    }
    
    template<class F3D>
    void nc_end_recv(int i, int j, int k, int X, int Y, int Z, F3D& F) const
    {
      float* p = static_cast<float*>(::end_recv_port(i,j,k, g_));
      if (p) {
	p++;                 /* Remote g->d##X */
	int face = (i+j+k) < 0 ? nx_[X] + 1 : 0;
	foreach_node(g_, X, face, [&](int x, int y, int z) { (&F(x,y,z).ex)[X] = *p++; });
      }
    }

    template<class F3D>
    void ec_end_recv(int i, int j, int k, int X, int Y, int Z, F3D& F) const
    {
      float* p = static_cast<float*>(::end_recv_port(i,j,k, g_));
      if (p) {
	p++;                 /* Remote g->d##X */
	int face = (i+j+k) < 0 ? nx_[X] + 1 : 0;
	foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Y] = *p++; });
	foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Z] = *p++; });
      }
    }

    //private:
    int nx_[3];
    const grid_t *g_;
  };

  void begin_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.ec_begin_recv(-1, 0, 0, 0,1,2);
    comm.ec_begin_recv( 0,-1, 0, 1,2,0);
    comm.ec_begin_recv( 0, 0,-1, 2,0,1);
    comm.ec_begin_recv( 1, 0, 0, 0,1,2);
    comm.ec_begin_recv( 0, 1, 0, 1,2,0);
    comm.ec_begin_recv( 0, 0, 1, 2,0,1);

    comm.ec_begin_send(-1, 0, 0, 0,1,2, F);
    comm.ec_begin_send( 0,-1, 0, 1,2,0, F);
    comm.ec_begin_send( 0, 0,-1, 2,0,1, F);
    comm.ec_begin_send( 1, 0, 0, 0,1,2, F);
    comm.ec_begin_send( 0, 1, 0, 1,2,0, F);
    comm.ec_begin_send( 0, 0, 1, 2,0,1, F);
  }

  void end_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);
    
    comm.ec_end_recv(-1, 0, 0, 0,1,2, F);
    comm.ec_end_recv( 0,-1, 0, 1,2,0, F);
    comm.ec_end_recv( 0, 0,-1, 2,0,1, F);
    comm.ec_end_recv( 1, 0, 0, 0,1,2, F);
    comm.ec_end_recv( 0, 1, 0, 1,2,0, F);
    comm.ec_end_recv( 0, 0, 1, 2,0,1, F);

    comm.end_send(-1, 0, 0);
    comm.end_send( 0,-1, 0);
    comm.end_send( 0, 0,-1);
    comm.end_send( 1, 0, 0);
    comm.end_send( 0, 1, 0);
    comm.end_send( 0, 0, 1);
  }

  void begin_remote_ghost_norm_e(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.nc_begin_recv(-1, 0, 0, 0,1,2);
    comm.nc_begin_recv( 0,-1, 0, 1,2,0);
    comm.nc_begin_recv( 0, 0,-1, 2,0,1);
    comm.nc_begin_recv( 1, 0, 0, 0,1,2);
    comm.nc_begin_recv( 0, 1, 0, 1,2,0);
    comm.nc_begin_recv( 0, 0, 1, 2,0,1);

    comm.nc_begin_send(-1, 0, 0, 0,1,2, F);
    comm.nc_begin_send( 0,-1, 0, 1,2,0, F);
    comm.nc_begin_send( 0, 0,-1, 2,0,1, F);
    comm.nc_begin_send( 1, 0, 0, 0,1,2, F);
    comm.nc_begin_send( 0, 1, 0, 1,2,0, F);
    comm.nc_begin_send( 0, 0, 1, 2,0,1, F);
  }

  void end_remote_ghost_norm_e(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.nc_end_recv(-1, 0, 0, 0,1,2, F);
    comm.nc_end_recv( 0,-1, 0, 1,2,0, F);
    comm.nc_end_recv( 0, 0,-1, 2,0,1, F);
    comm.nc_end_recv( 1, 0, 0, 0,1,2, F);
    comm.nc_end_recv( 0, 1, 0, 1,2,0, F);
    comm.nc_end_recv( 0, 0, 1, 2,0,1, F);
    
    comm.end_send(-1, 0, 0);
    comm.end_send( 0,-1, 0);
    comm.end_send( 0, 0,-1);
    comm.end_send( 1, 0, 0);
    comm.end_send( 0, 1, 0);
    comm.end_send( 0, 0, 1);
  }

  void begin_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int size, face, x, y, z;
    float *p;

# define BEGIN_RECV(i,j,k,X,Y,Z)				\
    begin_recv_port(i,j,k,(1+n##Y*n##Z)*sizeof(float),g)
    BEGIN_RECV(-1, 0, 0,x,y,z);
    BEGIN_RECV( 0,-1, 0,y,z,x);
    BEGIN_RECV( 0, 0,-1,z,x,y);
    BEGIN_RECV( 1, 0, 0,x,y,z);
    BEGIN_RECV( 0, 1, 0,y,z,x);
    BEGIN_RECV( 0, 0, 1,z,x,y);
# undef BEGIN_RECV

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {		\
      size = ( 1 + n##Y*n##Z )*sizeof(float);			\
      p = (float *)size_send_port( i, j, k, size, g );		\
      if( p ) {							\
	(*(p++)) = g->d##X;					\
	face = (i+j+k)<0 ? 1 : n##X;				\
	X##_FACE_LOOP(face) (*(p++)) = F(x,y,z).div_b_err;	\
	begin_send_port( i, j, k, size, g );			\
      }								\
    } END_PRIMITIVE
    BEGIN_SEND(-1, 0, 0,x,y,z);
    BEGIN_SEND( 0,-1, 0,y,z,x);
    BEGIN_SEND( 0, 0,-1,z,x,y);
    BEGIN_SEND( 1, 0, 0,x,y,z);
    BEGIN_SEND( 0, 1, 0,y,z,x);
    BEGIN_SEND( 0, 0, 1,z,x,y);
# undef BEGIN_SEND
  }

  void end_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int face, x, y, z;
    float *p, lw, rw;

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {		\
      p = (float *)end_recv_port(i,j,k,g);			\
      if( p ) {							\
	lw = (*(p++));                 /* Remote g->d##X */	\
	rw = (2.*g->d##X)/(lw+g->d##X);				\
	lw = (lw-g->d##X)/(lw+g->d##X);				\
	face = (i+j+k)<0 ? n##X+1 : 0; /* Interpolate */	\
	X##_FACE_LOOP(face)					\
	  F(x,y,z).div_b_err = rw*(*(p++)) +			\
	  lw*F(x+i,y+j,z+k).div_b_err;				\
      }								\
    } END_PRIMITIVE
    END_RECV(-1, 0, 0,x,y,z);
    END_RECV( 0,-1, 0,y,z,x);
    END_RECV( 0, 0,-1,z,x,y);
    END_RECV( 1, 0, 0,x,y,z);
    END_RECV( 0, 1, 0,y,z,x);
    END_RECV( 0, 0, 1,z,x,y);
# undef END_RECV

# define END_SEND(i,j,k,X,Y,Z) end_send_port(i,j,k,g)
    END_SEND(-1, 0, 0,x,y,z);
    END_SEND( 0,-1, 0,y,z,x);
    END_SEND( 0, 0,-1,z,x,y);
    END_SEND( 1, 0, 0,x,y,z);
    END_SEND( 0, 1, 0,y,z,x);
    END_SEND( 0, 0, 1,z,x,y);
# undef END_SEND
  }

};


#endif

