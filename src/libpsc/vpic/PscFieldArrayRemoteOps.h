
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

  void begin_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int size, face, x, y, z;
    float *p;

# define BEGIN_RECV(i,j,k,X,Y,Z)					\
    begin_recv_port(i,j,k,(1+n##Y*(n##Z+1)+n##Z*(n##Y+1))*sizeof(float),g)
    BEGIN_RECV(-1, 0, 0,x,y,z);
    BEGIN_RECV( 0,-1, 0,y,z,x);
    BEGIN_RECV( 0, 0,-1,z,x,y);
    BEGIN_RECV( 1, 0, 0,x,y,z);
    BEGIN_RECV( 0, 1, 0,y,z,x);
    BEGIN_RECV( 0, 0, 1,z,x,y);
# undef BEGIN_RECV

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {		\
      size = (1+n##Y*(n##Z+1)+n##Z*(n##Y+1))*sizeof(float);	\
      p = (float *)size_send_port( i, j, k, size, g );		\
      if( p ) {							\
	(*(p++)) = g->d##X;					\
	face = (i+j+k)<0 ? 1 : n##X;				\
	Z##Y##_EDGE_LOOP(face) (*(p++)) = F(x,y,z).cb##Y;	\
	Y##Z##_EDGE_LOOP(face) (*(p++)) = F(x,y,z).cb##Z;	\
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

  void end_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int face, x, y, z;
    float *p, lw, rw;

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {                        \
      p = (float *)end_recv_port(i,j,k,g);				\
      if( p ) {								\
	lw = (*(p++));                 /* Remote g->d##X */		\
	rw = (2.*g->d##X)/(lw+g->d##X);					\
	lw = (lw-g->d##X)/(lw+g->d##X);					\
	face = (i+j+k)<0 ? n##X+1 : 0; /* Interpolate */		\
	Z##Y##_EDGE_LOOP(face)						\
	  F(x,y,z).cb##Y = rw*(*(p++)) + lw*F(x+i,y+j,z+k).cb##Y;	\
	Y##Z##_EDGE_LOOP(face)						\
	  F(x,y,z).cb##Z = rw*(*(p++)) + lw*F(x+i,y+j,z+k).cb##Z;	\
      }									\
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

  void begin_remote_ghost_norm_e(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t *g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int size, face, x, y, z;
    float *p;

# define BEGIN_RECV(i,j,k,X,Y,Z)					\
    begin_recv_port(i,j,k,( 1 + (n##Y+1)*(n##Z+1) )*sizeof(float),g)
    BEGIN_RECV(-1, 0, 0,x,y,z);
    BEGIN_RECV( 0,-1, 0,y,z,x);
    BEGIN_RECV( 0, 0,-1,z,x,y);
    BEGIN_RECV( 1, 0, 0,x,y,z);
    BEGIN_RECV( 0, 1, 0,y,z,x);
    BEGIN_RECV( 0, 0, 1,z,x,y);
# undef BEGIN_RECV

# define BEGIN_SEND(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {		\
      size = ( 1+ (n##Y+1)*(n##Z+1) )*sizeof(float);		\
      p = (float *)size_send_port( i, j, k, size, g );		\
      if( p ) {							\
	(*(p++)) = g->d##X;					\
	face = (i+j+k)<0 ? 1 : n##X;				\
	X##_NODE_LOOP(face) (*(p++)) = F(x,y,z).e##X;	\
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

  void end_remote_ghost_norm_e(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    const grid_t* g = fa.g;
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    int face, x, y, z;
    float *p, lw, rw;

# define END_RECV(i,j,k,X,Y,Z) BEGIN_PRIMITIVE {			\
      p = (float *)end_recv_port(i,j,k,g);                              \
      if( p ) {                                                         \
	lw = (*(p++));                 /* Remote g->d##X */             \
	rw = (2.*g->d##X)/(lw+g->d##X);                                 \
	lw = (lw-g->d##X)/(lw+g->d##X);                                 \
	face = (i+j+k)<0 ? n##X+1 : 0; /* Interpolate */                \
	X##_NODE_LOOP(face)                                             \
	  F(x,y,z).e##X = rw*(*(p++)) + lw*F(x+i,y+j,z+k).e##X;		\
      }                                                                 \
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

