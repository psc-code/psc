
#ifndef PSC_PARTICLES_OPS
#define PSC_PARTICLES_OPS

#define HAS_V4_PIPELINE

template<class P>
struct PscParticlesOps {
  typedef P Particles;
  typedef typename Particles::Species Species;
  typedef typename Particles::FieldArray FieldArray;
  typedef typename Particles::Interpolator Interpolator;
  typedef typename Particles::Accumulator Accumulator;
  typedef typename Accumulator::Block AccumulatorBlock;
  
  PscParticlesOps(vpic_simulation *simulation) : simulation_(simulation) { }

  // ----------------------------------------------------------------------
  // move_p
  //
  // move_p moves the particle m->p by m->dispx, m->dispy, m->dispz
  // depositing particle current as it goes. If the particle was moved
  // sucessfully (particle mover is no longer in use) returns 0. If the
  // particle interacted with something this routine could not handle,
  // this routine returns 1 (particle mover is still in use). On a
  // successful move, the particle position is updated and m->dispx,
  // m->dispy and m->dispz are zerod. On a partial move, the particle
  // position is updated to the point where the particle interacted and
  // m->dispx, m->dispy, m->dispz contains the remaining particle
  // displacement. The displacements are the physical displacments
  // normalized current cell size.
  //
  // Because move_p is frequently called, it does not check its input
  // arguments. Higher level routines are responsible for insuring valid
  // arguments.
  //
  // Note: changes here likely need to be reflected in SPE accelerated
  // version as well.

#define ACCUMULATE_J(u,d,X,Y,Z,offset)					\
      v4  = q*u##X;   /* v2 = q ux                            */	\
      v1  = v4*d##Y;  /* v1 = q ux dy                         */	\
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */	\
      v1 += v4;       /* v1 = q ux (1+dy)                     */	\
      v4  = one+d##Z; /* v4 = 1+dz                            */	\
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */	\
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */	\
      v4  = one-d##Z; /* v4 = 1-dz                            */	\
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */	\
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */	\
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */	\
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */	\
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */	\
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */	\
      a[offset+0] += v0;						\
      a[offset+1] += v1;						\
      a[offset+2] += v2;						\
      a[offset+3] += v3

#if defined(V4_ACCELERATION)

  // High performance variant based on SPE accelerated version

  int move_p(particle_t       * RESTRICT ALIGNED(128) p,
	     particle_mover_t * RESTRICT ALIGNED(16)  pm,
	     AccumulatorBlock acc_block,
	     const grid_t     *                       g,
	     const float                              qsp)
  {

    using namespace v4;

    /*const*/ v4float one( 1.f );
    /*const*/ v4float tiny( 1e-37f );
    /*const*/ v4int   sign_bits( 1<<31 );

    v4float dr, r, u, q, q3;
    v4float sgn_dr, s, sdr;
    v4float v0, v1, v2, v3, v4, v5, _stack_vf;
    v4int bits, _stack_vi;

    float * RESTRICT ALIGNED(16) stack_vf = (float *)&_stack_vf;
    int   * RESTRICT ALIGNED(16) stack_vi =   (int *)&_stack_vi;
    float f0, f1;
    int32_t n, voxel;
    int64_t neighbor;
    int type;

    load_4x1( &pm->dispx, dr );  n     = pm->i;
    load_4x1( &p[n].dx,   r  );  voxel = p[n].i;
    load_4x1( &p[n].ux,   u  );

    q  = v4float(qsp)*splat<3>(u); // q  = p_q,   p_q,   p_q,   D/C
    q3 = v4float(1.f/3.f)*q;      // q3 = p_q/3, p_q/3, p_q/3, D/C
    dr = shuffle<0,1,2,2>( dr );  // dr = p_ddx, p_ddy, p_ddz, D/C 
    r  = shuffle<0,1,2,2>( r );  // r  = p_dx,  p_dy,  p_dz,  D/C
  
    for(;;) {

      // At this point:
      //   r     = current particle position in local voxel coordinates
      //           (note the current voxel is on [-1,1]^3.
      //   dr    = remaining particle displacment
      //           (note: this is in voxel edge lengths!)
      //   voxel = local voxel of particle
      // Thus, in the local coordinate system, it is desired to move the
      // particle through all points in the local coordinate system:
      //   streak_r(s) = r + 2 disp s for s in [0,1]
      //
      // Determine the fractional length and type of current
      // streak made by the particle through this voxel.  The streak
      // ends on either the first voxel face intersected by the
      // particle track or at the end of the particle track.
      //
      // Note: a divide by zero cannot occur below due to the shift of
      // the denominator by tiny.  Also, the shift by tiny is large
      // enough that the divide will never overflow when dr is tiny
      // (|sgn_dr-r|<=2 => 2/tiny = 2e+37 < FLT_MAX = 3.4e38).
      // Likewise, due to speed of light limitations, generally dr
      // cannot get much larger than 1 or so and the numerator, if not
      // zero, can generally never be smaller than FLT_EPS/2.  Thus,
      // likewise, the divide will never underflow either. 

      // FIXME: THIS COULD PROBABLY BE DONE EVEN FASTER 
      sgn_dr = copysign( one,  dr );
      v0     = copysign( tiny, dr );
      store_4x1( (sgn_dr-r) / ((dr+dr)+v0), stack_vf );
      /**/                          type = 3;             f0 = 1;
      f1 = stack_vf[0]; if( f1<f0 ) type = 0; if( f1<f0 ) f0 = f1; // Branchless cmov 
      f1 = stack_vf[1]; if( f1<f0 ) type = 1; if( f1<f0 ) f0 = f1;
      f1 = stack_vf[2]; if( f1<f0 ) type = 2; if( f1<f0 ) f0 = f1;
      s = v4float( f0 );

      // At this point:
      //   type = 0,1 or 2 ... streak ends on a x,y or z-face respectively
      //          3        ... streak ends at end of the particle track
      //   s    = SPLAT( normalized length of the current streak )
      //   sgn_dr indicates the sign streak displacements.  This is
      //     useful to determine whether or not the streak hit a face.
      //
      // Compute the streak midpoint the normalized displacement,
      // update the particle position and remaining displacment,
      // compute accumulator coefficients, finish up fetching the
      // voxel needed for this streak and accumule the streak.  Note:
      // accumulator values are 4 times the total physical charge that
      // passed through the appropriate current quadrant in a
      // timestep.

      sdr = s*dr;      // sdr = ux,          uy,          uz,          D/C
      v5  = r + sdr;   // v5  = dx,          dy,          dz,          D/C

      dr -= sdr;       // dr  = p_ddx',      p_ddy',      p_ddz',      D/C
      r  += sdr + sdr; // r   = p_dx',       p_dy',       p_dz',       D/C

      v4  = q*sdr;     // v4  = q ux,        q uy,        q uz,        D/C
      v1  = v4*shuffle<1,2,0,3>( v5 );
      /**/             // v1  = q ux dy,     q uy dz,     q uz dx,     D/C
      v0  = v4 - v1;   // v0  = q ux(1-dy),  q uy(1-dz),  q uz(1-dx),  D/C
      v1 += v4;        // v1  = q ux(1+dy),  q uy(1+dz),  q uz(1+dx),  D/C

      v5  = shuffle<2,0,1,3>( v5 ); // v5 = dz, dx, dy, D/C
      v4  = one + v5;  // v4  = 1+dz,        1+dx,        1+dy,        D/C
      v2  = v0*v4;     // v2  = q ux(1-dy)(1+dz), ...,                 D/C
      v3  = v1*v4;     // v3  = q ux(1+dy)(1+dz), ...,                 D/C
      v4  = one - v5;  // v4  = 1-dz,        1-dx,        1-dy,        D/C
      v0 *= v4;        // v0  = q ux(1-dy)(1-dz), ...,                 D/C
      v1 *= v4;        // v1  = q ux(1+dy)(1-dz), ...,                 D/C

      //v4  = ((q3*splat(sdr,0))*splat(sdr,1))*splat(sdr,2);
      v4  = ((q3*splat<0>(sdr))*splat<1>(sdr))*splat<2>(sdr);
      // FIXME: splat ambiguity in v4 prevents flattening
      /**/             // v4  = q ux uy uz/3,q ux uy uz/3,q ux uy uz/3,D/C
      v0 += v4;        // v0  = q ux[(1-dy)(1-dz)+uy uz/3], ...,       D/C
      v1 -= v4;        // v1  = q ux[(1+dy)(1-dz)-uy uz/3], ...,       D/C
      v2 -= v4;        // v2  = q ux[(1-dy)(1+dz)-uy uz/3], ...,       D/C
      v3 += v4;        // v3  = q ux[(1+dy)(1+dz)+uy uz/3], ...,       D/C

      transpose( v0, v1, v2, v3 );

      increment_4x1( acc_block[voxel].jx, v0 );
      increment_4x1( acc_block[voxel].jy, v1 );
      increment_4x1( acc_block[voxel].jz, v2 );

      // If streak ended at the end of the particle track, this mover
      // was succesfully processed.  Should be just under ~50% of the
      // time.
       
      if( type==3 ) { store_4x1( r, &p[n].dx ); p[n].i = voxel; break; }

      // Streak terminated on a voxel face.  Determine if the particle
      // crossed into a local voxel or if it hit a boundary.  Convert
      // the particle coordinates accordingly.  Note: Crossing into a
      // local voxel should happen the most of the time once we get to
      // this point; hitting a structure or parallel domain boundary
      // should usually be a rare event. */

      clear_4x1( stack_vi );
      stack_vi[type] = -1;
      load_4x1( stack_vi, bits );
      r = merge( bits, sgn_dr, r ); // Avoid roundoff fiascos--put the
      // particle _exactly_ on the
      // boundary.
      bits &= sign_bits; // bits(type)==(-0.f) and 0 elsewhere

      // Determine if the particle crossed into a local voxel or if it
      // hit a boundary.  Convert the particle coordinates accordingly.
      // Note: Crossing into a local voxel should happen the other ~50%
      // of time; hitting a structure and parallel domain boundary
      // should usually be a rare event.  Note: the entry / exit
      // coordinate for the particle is guaranteed to be +/-1 _exactly_
      // for the particle.

      store_4x1( sgn_dr, stack_vf ); if( stack_vf[type]>0 ) type += 3;
      neighbor = g->neighbor[ 6*voxel + type ];

      if( UNLIKELY( neighbor==reflect_particles ) ) {

	// Hit a reflecting boundary condition.  Reflect the particle
	// momentum and remaining displacement and keep moving the
	// particle.

	dr = toggle_bits( bits, dr );
	u  = toggle_bits( bits, u  );
	store_4x1( u, &p[n].ux );
	continue;
      }

      if( UNLIKELY( neighbor<g->rangel || neighbor>g->rangeh ) ) {

	// Cannot handle the boundary condition here.  Save the updated
	// particle position and update the remaining displacement in
	// the particle mover.

	store_4x1( r, &p[n].dx );    p[n].i = 8*voxel + type;
	store_4x1( dr, &pm->dispx ); pm->i  = n;
	return 1; // Mover still in use
      }

      // Crossed into a normal voxel.  Update the voxel index, convert the
      // particle coordinate system and keep moving the particle.
      
      voxel = (int32_t)( neighbor - g->rangel );
      r = toggle_bits( bits, r );
    }

    return 0; // Mover not in use
  }

#else

  int move_p(particle_t       * ALIGNED(128) p0,
	     particle_mover_t * ALIGNED(16)  pm,
	     AccumulatorBlock acc_block,
	     const grid_t     *              g,
	     const float                     qsp )
  {
    float s_midx, s_midy, s_midz;
    float s_dispx, s_dispy, s_dispz;
    float s_dir[3];
    float v0, v1, v2, v3, v4, v5, q;
    int axis, face;
    int64_t neighbor;
    float *a;
    particle_t * ALIGNED(32) p = p0 + pm->i;

    // FIXME, this uses double precision constants in a bunch of places
    const float one = 1.;
    
    q = qsp*p->w;

    for(;;) {
      s_midx = p->dx;
      s_midy = p->dy;
      s_midz = p->dz;

      s_dispx = pm->dispx;
      s_dispy = pm->dispy;
      s_dispz = pm->dispz;

      s_dir[0] = (s_dispx>0) ? 1 : -1;
      s_dir[1] = (s_dispy>0) ? 1 : -1;
      s_dir[2] = (s_dispz>0) ? 1 : -1;
    
      // Compute the twice the fractional distance to each potential
      // streak/cell face intersection.
      v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx;
      v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
      v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

      // Determine the fractional length and axis of current streak. The
      // streak ends on either the first face intersected by the
      // particle track or at the end of the particle track.
      // 
      //   axis 0,1 or 2 ... streak ends on a x,y or z-face respectively
      //   axis 3        ... streak ends at end of the particle track
      /**/      v3=2,  axis=3;
      if(v0<v3) v3=v0, axis=0;
      if(v1<v3) v3=v1, axis=1;
      if(v2<v3) v3=v2, axis=2;
      v3 *= 0.5;

      // Compute the midpoint and the normalized displacement of the streak
      s_dispx *= v3;
      s_dispy *= v3;
      s_dispz *= v3;
      s_midx += s_dispx;
      s_midy += s_dispy;
      s_midz += s_dispz;

      // Accumulate the streak.  Note: accumulator values are 4 times
      // the total physical charge that passed through the appropriate
      // current quadrant in a time-step
      v5 = q*s_dispx*s_dispy*s_dispz*(1./3.);
      a = (float *) &acc_block[p->i];
      
      ACCUMULATE_J(s_disp,s_mid,x,y,z, 0);
      ACCUMULATE_J(s_disp,s_mid,y,z,x, 4);
      ACCUMULATE_J(s_disp,s_mid,z,x,y, 8);

      // Compute the remaining particle displacment
      pm->dispx -= s_dispx;
      pm->dispy -= s_dispy;
      pm->dispz -= s_dispz;

      // Compute the new particle offset
      p->dx += s_dispx+s_dispx;
      p->dy += s_dispy+s_dispy;
      p->dz += s_dispz+s_dispz;

      // If an end streak, return success (should be ~50% of the time)

      if( axis==3 ) break;

      // Determine if the particle crossed into a local cell or if it
      // hit a boundary and convert the coordinate system accordingly.
      // Note: Crossing into a local cell should happen ~50% of the
      // time; hitting a boundary is usually a rare event.  Note: the
      // entry / exit coordinate for the particle is guaranteed to be
      // +/-1 _exactly_ for the particle.

      v0 = s_dir[axis];
      (&(p->dx))[axis] = v0; // Avoid roundoff fiascos--put the particle
      // _exactly_ on the boundary.
      face = axis; if( v0>0 ) face += 3;
      neighbor = g->neighbor[ 6*p->i + face ];
    
      if( UNLIKELY( neighbor==reflect_particles ) ) {
	// Hit a reflecting boundary condition.  Reflect the particle
	// momentum and remaining displacement and keep moving the
	// particle.
	(&(p->ux    ))[axis] = -(&(p->ux    ))[axis];
	(&(pm->dispx))[axis] = -(&(pm->dispx))[axis];
	continue;
      }

      if( UNLIKELY( neighbor<g->rangel || neighbor>g->rangeh ) ) {
	// Cannot handle the boundary condition here.  Save the updated
	// particle position, face it hit and update the remaining
	// displacement in the particle mover.
	p->i = 8*p->i + face;
	return 1; // Return "mover still in use"
      }

      // Crossed into a normal voxel.  Update the voxel index, convert the
      // particle coordinate system and keep moving the particle.
    
      p->i = neighbor - g->rangel; // Compute local index of neighbor
      /**/                         // Note: neighbor - g->rangel < 2^31 / 6
      (&(p->dx))[axis] = -v0;      // Convert coordinate system
    }

    return 0; // Return "mover not in use"
  }

#endif

  // ----------------------------------------------------------------------
  // inject_particles
  
  void inject_particle(Particles& vmprts, Accumulator& accumulator, FieldArray& fa,
		       const struct psc_particle_inject *prt)
  {
    auto sp = vmprts.find(prt->kind);

    double x = prt->x[0], y = prt->x[1], z = prt->x[2];
    double ux = prt->u[0], uy = prt->u[1], uz = prt->u[2];
    double w = prt->w, age = 0.;
    int update_rhob = 0;

    int ix, iy, iz;

    grid_t *grid = sp->g;
    const double x0 = (double)grid->x0, y0 = (double)grid->y0, z0 = (double)grid->z0;
    const double x1 = (double)grid->x1, y1 = (double)grid->y1, z1 = (double)grid->z1;
    const int    nx = grid->nx,         ny = grid->ny,         nz = grid->nz;
    
    // Do not inject if the particle is strictly outside the local domain
    // or if a far wall of local domain shared with a neighbor
    
    if ((x<x0) | (x>x1) | ( (x==x1) & (grid->bc[BOUNDARY(1,0,0)]>=0))) return;
    if ((y<y0) | (y>y1) | ( (y==y1) & (grid->bc[BOUNDARY(0,1,0)]>=0))) return;
    if ((z<z0) | (z>z1) | ( (z==z1) & (grid->bc[BOUNDARY(0,0,1)]>=0))) return;
    
    // This node should inject the particle
    
    if (sp->np>=sp->max_np) ERROR(( "No room to inject particle" ));
    
    // Compute the injection cell and coordinate in cell coordinate system
    // BJA:  Note the use of double precision here for accurate particle 
    //       placement on large meshes. 
 
    // The ifs allow for injection on the far walls of the local computational
    // domain when necessary
    
    x  = ((double)nx)*((x-x0)/(x1-x0)); // x is rigorously on [0,nx]
    ix = (int)x;                        // ix is rigorously on [0,nx]
    x -= (double)ix;                    // x is rigorously on [0,1)
    x  = (x+x)-1;                       // x is rigorously on [-1,1)
    if( ix==nx ) x = 1;                 // On far wall ... conditional move
    if( ix==nx ) ix = nx-1;             // On far wall ... conditional move
    ix++;                               // Adjust for mesh indexing

    y  = ((double)ny)*((y-y0)/(y1-y0)); // y is rigorously on [0,ny]
    iy = (int)y;                        // iy is rigorously on [0,ny]
    y -= (double)iy;                    // y is rigorously on [0,1)
    y  = (y+y)-1;                       // y is rigorously on [-1,1)
    if( iy==ny ) y = 1;                 // On far wall ... conditional move
    if( iy==ny ) iy = ny-1;             // On far wall ... conditional move
    iy++;                               // Adjust for mesh indexing

    z  = ((double)nz)*((z-z0)/(z1-z0)); // z is rigorously on [0,nz]
    iz = (int)z;                        // iz is rigorously on [0,nz]
    z -= (double)iz;                    // z is rigorously on [0,1)
    z  = (z+z)-1;                       // z is rigorously on [-1,1)
    if( iz==nz ) z = 1;                 // On far wall ... conditional move
    if( iz==nz ) iz = nz-1;             // On far wall ... conditional move
    iz++;                               // Adjust for mesh indexing

    particle_t * p = sp->p + (sp->np++);
    p->dx = (float)x; // Note: Might be rounded to be on [-1,1]
    p->dy = (float)y; // Note: Might be rounded to be on [-1,1]
    p->dz = (float)z; // Note: Might be rounded to be on [-1,1]
    p->i  = VOXEL(ix,iy,iz, nx,ny,nz);
    p->ux = (float)ux;
    p->uy = (float)uy;
    p->uz = (float)uz;
    p->w  = w;

    if (update_rhob) accumulate_rhob(fa, p, -sp->q);

    if (age!=0) {
      if( sp->nm >= sp->max_nm )
	WARNING(( "No movers available to age injected  particle" ));
      particle_mover_t * pm = sp->pm + sp->nm;
      age *= grid->cvac*grid->dt/sqrt( ux*ux + uy*uy + uz*uz + 1 );
      pm->dispx = ux*age*grid->rdx;
      pm->dispy = uy*age*grid->rdy;
      pm->dispz = uz*age*grid->rdz;
      pm->i     = sp->np-1;
      sp->nm += move_p( sp->p, pm, accumulator[0], grid, sp->q );
    }
    
  }
  
  // ----------------------------------------------------------------------
  // advance_p

  typedef struct particle_mover_seg {
    int nm;                             // Number of movers used
    int n_ignored;                      // Number of movers ignored
  } particle_mover_seg_t;

  void advance_p_pipeline(typename Particles::iterator sp,
			  AccumulatorBlock acc_block,
			  Interpolator& interpolator,
			  particle_mover_seg_t *seg,
			  particle_t * ALIGNED(128) p, int n,
			  particle_mover_t * ALIGNED(16) pm, int max_nm)
  {
    particle_t           * ALIGNED(128) p0 = sp->p;
    const grid_t *                      g  = sp->g;

    const interpolator_t * ALIGNED(16)  f;
    float                * ALIGNED(16)  a;

    const float qdt_2mc  = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
    const float cdt_dx   = sp->g->cvac*sp->g->dt*sp->g->rdx;
    const float cdt_dy   = sp->g->cvac*sp->g->dt*sp->g->rdy;
    const float cdt_dz   = sp->g->cvac*sp->g->dt*sp->g->rdz;
    const float qsp      = sp->q;

    const float one            = 1.;
    const float one_third      = 1./3.;
    const float two_fifteenths = 2./15.;

    float dx, dy, dz, ux, uy, uz, q;
    float hax, hay, haz, cbx, cby, cbz;
    float v0, v1, v2, v3, v4, v5;

    int ii;
  
    DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

    int nm = 0;
    int n_ignored = 0;

    // Process particles for this pipeline

    for(;n;n--,p++) {
      dx   = p->dx;                             // Load position
      dy   = p->dy;
      dz   = p->dz;
      ii   = p->i;
      f    = &interpolator[ii];                 // Interpolate E
      hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
			  dz*( f->dexdz + dy*f->d2exdydz ) );
      hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
			  dx*( f->deydx + dz*f->d2eydzdx ) );
      haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
			  dy*( f->dezdy + dx*f->d2ezdxdy ) );
      cbx  = f->cbx + dx*f->dcbxdx;             // Interpolate B
      cby  = f->cby + dy*f->dcbydy;
      cbz  = f->cbz + dz*f->dcbzdz;
      ux   = p->ux;                             // Load momentum
      uy   = p->uy;
      uz   = p->uz;
      q    = p->w;
      ux  += hax;                               // Half advance E
      uy  += hay;
      uz  += haz;
      v0   = qdt_2mc/sqrtf(one + (ux*ux + (uy*uy + uz*uz)));
      /**/                                      // Boris - scalars
      v1   = cbx*cbx + (cby*cby + cbz*cbz);
      v2   = (v0*v0)*v1;
      v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
      v4   = v3/(one+v1*(v3*v3));
      v4  += v4;
      v0   = ux + v3*( uy*cbz - uz*cby );       // Boris - uprime
      v1   = uy + v3*( uz*cbx - ux*cbz );
      v2   = uz + v3*( ux*cby - uy*cbx );
      ux  += v4*( v1*cbz - v2*cby );            // Boris - rotation
      uy  += v4*( v2*cbx - v0*cbz );
      uz  += v4*( v0*cby - v1*cbx );
      ux  += hax;                               // Half advance E
      uy  += hay;
      uz  += haz;
      p->ux = ux;                               // Store momentum
      p->uy = uy;
      p->uz = uz;
      v0   = one/sqrtf(one + (ux*ux+ (uy*uy + uz*uz)));
      /**/                                      // Get norm displacement
      ux  *= cdt_dx;
      uy  *= cdt_dy;
      uz  *= cdt_dz;
      ux  *= v0;
      uy  *= v0;
      uz  *= v0;
      v0   = dx + ux;                           // Streak midpoint (inbnds)
      v1   = dy + uy;
      v2   = dz + uz;
      v3   = v0 + ux;                           // New position
      v4   = v1 + uy;
      v5   = v2 + uz;

      // FIXME-KJB: COULD SHORT CIRCUIT ACCUMULATION IN THE CASE WHERE QSP==0!
      if(  v3<=one &&  v4<=one &&  v5<=one &&   // Check if inbnds
	   -v3<=one && -v4<=one && -v5<=one ) {

	// Common case (inbnds).  Note: accumulator values are 4 times
	// the total physical charge that passed through the appropriate
	// current quadrant in a time-step

	q *= qsp;
	p->dx = v3;                             // Store new position
	p->dy = v4;
	p->dz = v5;
	dx = v0;                                // Streak midpoint
	dy = v1;
	dz = v2;
	v5 = q*ux*uy*uz*one_third;              // Compute correction
	a  = (float *) &acc_block[ii];          // Get accumulator

	ACCUMULATE_J(u,d, x,y,z, 0);
	ACCUMULATE_J(u,d, y,z,x, 4);
	ACCUMULATE_J(u,d, z,x,y, 8);

#     undef ACCUMULATE_J

      } else {                                    // Unlikely
	local_pm->dispx = ux;
	local_pm->dispy = uy;
	local_pm->dispz = uz;
	local_pm->i     = p - p0;

	if( move_p( p0, local_pm, acc_block, g, qsp ) ) { // Unlikely
	  if( nm<max_nm ) {
	    pm[nm++] = local_pm[0];
	  }
	  else {
	    n_ignored++;                 // Unlikely
	  } // if
	} // if
      }

    }

    seg->nm        = nm;
    seg->n_ignored = n_ignored;
  }

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  void advance_p_pipeline_v4(typename Particles::iterator sp,
			     AccumulatorBlock acc_block,
			     Interpolator& interpolator,
			     particle_mover_seg_t *seg,
			     particle_t * ALIGNED(128) p, int n,
			     particle_mover_t * ALIGNED(16) pm, int max_nm)
  {
    using namespace v4;

    particle_t           * ALIGNED(128) p0 = sp->p;
    const grid_t *                      g  = sp->g;

    float                * ALIGNED(16)  vp0;
    float                * ALIGNED(16)  vp1;
    float                * ALIGNED(16)  vp2;
    float                * ALIGNED(16)  vp3;

    const v4float qdt_2mc  = (sp->q*sp->g->dt)/(2*sp->m*sp->g->cvac);
    const v4float cdt_dx   = sp->g->cvac*sp->g->dt*sp->g->rdx;
    const v4float cdt_dy   = sp->g->cvac*sp->g->dt*sp->g->rdy;
    const v4float cdt_dz   = sp->g->cvac*sp->g->dt*sp->g->rdz;
    const v4float qsp      = sp->q;

    const v4float one = 1.;
    const v4float one_third = 1./3.;
    const v4float two_fifteenths = 2./15.;
    const v4float neg_one = -1.;

    const float _qsp = sp->q;

    v4float dx, dy, dz, ux, uy, uz, q;
    v4float hax, hay, haz, cbx, cby, cbz;
    v4float v0, v1, v2, v3, v4, v5;
    v4int   ii, outbnd;

    DECLARE_ALIGNED_ARRAY( particle_mover_t, 16, local_pm, 1 );

    int nq = n >> 2;
  
    int nm   = 0;
    int n_ignored = 0;

    // Process the particle quads for this pipeline

    for( ; nq; nq--, p+=4 ) {
      load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

      // Interpolate fields
      vp0 = (float * ALIGNED(16)) (&interpolator[ii(0)]);
      vp1 = (float * ALIGNED(16)) (&interpolator[ii(1)]);
      vp2 = (float * ALIGNED(16)) (&interpolator[ii(2)]);
      vp3 = (float * ALIGNED(16)) (&interpolator[ii(3)]);
      load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( v2, dy, v1 ), dz, fma( v0, dy, hax ) );
      load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( v5, dz, v4 ), dx, fma( v3, dz, hay ) );
      load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( v2, dx, v1 ), dy, fma( v0, dx, haz ) );
      load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
      /**/                                                    cby = fma( v4, dy, cby );
      load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

      // Update momentum
      // If you are willing to eat a 5-10% performance hit,
      // v0 = qdt_2mc/sqrt(blah) is a few ulps more accurate (but still
      // quite in the noise numerically) for cyclotron frequencies
      // approaching the nyquist frequency.

      load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
      ux += hax;
      uy += hay;
      uz += haz;
      v0  = qdt_2mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
      v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
      v2  = (v0*v0)*v1;
      v3  = v0*fma( fma( two_fifteenths, v2, one_third ), v2, one );
      v4  = v3*rcp(fma( v3*v3, v1, one ));
      v4 += v4;
      v0  = fma( fms( uy,cbz, uz*cby ), v3, ux );
      v1  = fma( fms( uz,cbx, ux*cbz ), v3, uy );
      v2  = fma( fms( ux,cby, uy*cbx ), v3, uz );
      ux  = fma( fms( v1,cbz, v2*cby ), v4, ux );
      uy  = fma( fms( v2,cbx, v0*cbz ), v4, uy );
      uz  = fma( fms( v0,cby, v1*cbx ), v4, uz );
      ux += hax;
      uy += hay;
      uz += haz;
      store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
    
      // Update the position of inbnd particles
      v0  = rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
      ux *= cdt_dx;
      uy *= cdt_dy;
      uz *= cdt_dz;
      ux *= v0;
      uy *= v0;
      uz *= v0;      // ux,uy,uz are normalized displ (relative to cell size)
      v0  = dx + ux;
      v1  = dy + uy;
      v2  = dz + uz; // New particle midpoint
      v3  = v0 + ux;
      v4  = v1 + uy;
      v5  = v2 + uz; // New particle position
      outbnd = (v3>one) | (v3<neg_one) |
	(v4>one) | (v4<neg_one) |
	(v5>one) | (v5<neg_one);
      v3  = merge(outbnd,dx,v3); // Do not update outbnd particles
      v4  = merge(outbnd,dy,v4);
      v5  = merge(outbnd,dz,v5);
      store_4x4_tr(v3,v4,v5,ii,&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx);
    
      // Accumulate current of inbnd particles
      // Note: accumulator values are 4 times the total physical charge that
      // passed through the appropriate current quadrant in a time-step
      q  = czero(outbnd,q*qsp);       // Do not accumulate outbnd particles
      dx = v0;                       // Streak midpoint (valid for inbnd only)
      dy = v1;
      dz = v2;
      v5 = q*ux*uy*uz*one_third;     // Charge conservation correction
      // seems to cause a measurable slow-down
      // FIXME, at least part of the solution should be to pass in a view of
      // the accumulator array for the current block, which also would mesh better
      // with only having a single block in the future at some point
      vp0 = (float * ALIGNED(16)) &acc_block[ii(0)]; // Accumlator pointers
      vp1 = (float * ALIGNED(16)) &acc_block[ii(1)];
      vp2 = (float * ALIGNED(16)) &acc_block[ii(2)];
      vp3 = (float * ALIGNED(16)) &acc_block[ii(3)];

#   define ACCUMULATE_J(X,Y,Z,offset)					\
      v4  = q*u##X;   /* v4 = q ux                            */	\
      v1  = v4*d##Y;  /* v1 = q ux dy                         */	\
      v0  = v4-v1;    /* v0 = q ux (1-dy)                     */	\
      v1 += v4;       /* v1 = q ux (1+dy)                     */	\
      v4  = one+d##Z; /* v4 = 1+dz                            */	\
      v2  = v0*v4;    /* v2 = q ux (1-dy)(1+dz)               */	\
      v3  = v1*v4;    /* v3 = q ux (1+dy)(1+dz)               */	\
      v4  = one-d##Z; /* v4 = 1-dz                            */	\
      v0 *= v4;       /* v0 = q ux (1-dy)(1-dz)               */	\
      v1 *= v4;       /* v1 = q ux (1+dy)(1-dz)               */	\
      v0 += v5;       /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */	\
      v1 -= v5;       /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */	\
      v2 -= v5;       /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */	\
      v3 += v5;       /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */	\
      transpose(v0,v1,v2,v3);						\
      increment_4x1(vp0+offset,v0);					\
      increment_4x1(vp1+offset,v1);					\
      increment_4x1(vp2+offset,v2);					\
      increment_4x1(vp3+offset,v3)

      ACCUMULATE_J( x,y,z, 0 );
      ACCUMULATE_J( y,z,x, 4 );
      ACCUMULATE_J( z,x,y, 8 );

#   undef ACCUMULATE_J

      // Update position and accumulate outbnd

#   define MOVE_OUTBND(N)                                               \
      if( outbnd(N) ) {                       /* Unlikely */		\
	local_pm->dispx = ux(N);					\
	local_pm->dispy = uy(N);					\
	local_pm->dispz = uz(N);					\
	local_pm->i     = (p - p0) + N;					\
	if( move_p( p0, local_pm, acc_block, g, _qsp ) ) { /* Unlikely */	\
      if( nm<max_nm ) copy_4x1( &pm[nm++], local_pm );			\
      else            n_ignored++;             /* Unlikely */		\
    }									\
    }

      MOVE_OUTBND(0);
      MOVE_OUTBND(1);
      MOVE_OUTBND(2);
      MOVE_OUTBND(3);

#   undef MOVE_OUTBND

    }

    seg->nm        = nm;
    seg->n_ignored = n_ignored;
  }

#endif
          
  void advance_p(typename Particles::iterator sp, Accumulator& accumulator, Interpolator& interpolator)
  {
    DECLARE_ALIGNED_ARRAY( particle_mover_seg_t, 128, seg, 1 );

    sp->nm = 0;

    particle_t *p = sp->p;
    int n = sp->np & ~15;
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)
    advance_p_pipeline_v4(sp, accumulator[1], interpolator, seg, p, n,
			  sp->pm + sp->nm, sp->max_nm - sp->nm);
#else
    advance_p_pipeline(sp, accumulator[1], interpolator, seg, p, n,
		       sp->pm + sp->nm, sp->max_nm - sp->nm);
#endif
    sp->nm += seg->nm;

    if (seg->n_ignored)
      WARNING(( "Pipeline %i ran out of storage for %i movers",
		0, seg->n_ignored ));
  
    p += n;
    n = sp->np - n;
    advance_p_pipeline(sp, accumulator[0], interpolator, seg, p, n,
		       sp->pm + sp->nm, sp->max_nm - sp->nm);
    sp->nm += seg->nm;

    if (seg->n_ignored)
      WARNING(( "Pipeline %i ran out of storage for %i movers",
		1, seg->n_ignored ));
  }


  void advance_p(Particles& vmprts, Accumulator& accumulator,
		 Interpolator& interpolator)
  {
    for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
      TIC advance_p(sp, accumulator, interpolator); TOC(advance_p, 1);
    }
  }

  // ----------------------------------------------------------------------
  // boundary_p

  void
  boundary_p_(particle_bc_t* pbc_list, Particles& vmprts, FieldArray& fa,
	      AccumulatorBlock acc_block)
  {
#ifdef V4_ACCELERATION
    using namespace v4;
#endif
    
    enum { MAX_PBC = 32, MAX_SP = 32 };

    // Gives the local mp port associated with a local face
    static const int f2b[6]  = { BOUNDARY(-1, 0, 0),
				 BOUNDARY( 0,-1, 0),
				 BOUNDARY( 0, 0,-1),
				 BOUNDARY( 1, 0, 0),
				 BOUNDARY( 0, 1, 0),
				 BOUNDARY( 0, 0, 1) };

    // Gives the remote mp port associated with a local face
    static const int f2rb[6] = { BOUNDARY( 1, 0, 0),
				 BOUNDARY( 0, 1, 0),
				 BOUNDARY( 0, 0, 1),
				 BOUNDARY(-1, 0, 0),
				 BOUNDARY( 0,-1, 0),
				 BOUNDARY( 0, 0,-1) };

    // Gives the axis associated with a local face
    static const int axis[6]  = { 0, 1, 2,  0,  1,  2 };

    // Gives the location of sending face on the receiver
    static const float dir[6] = { 1, 1, 1, -1, -1, -1 };

    // Temporary store for local particle injectors
    // FIXME: Ugly static usage
    static particle_injector_t * RESTRICT ALIGNED(16) ci = NULL;
    static int max_ci = 0;

    int n_send[6], n_recv[6], n_ci;

    int face;

    // Check input args 

    if (vmprts.empty()) return; // Nothing to do if no species 

    assert(num_particle_bc(pbc_list) == 0);

    // Unpack fields

    grid_t  * RESTRICT              g = fa.g;

    // Unpack the grid

    const int64_t * RESTRICT ALIGNED(128) neighbor = g->neighbor;
    /**/  mp_t    * RESTRICT              mp       = g->mp;
    const int64_t rangel = g->rangel;
    const int64_t rangeh = g->rangeh;
    const int64_t rangem = g->range[world_size];
    /*const*/ int bc[6], shared[6];
    /*const*/ int64_t range[6];
    for( face=0; face<6; face++ ) {
      bc[face] = g->bc[f2b[face]];
      shared[face] = (bc[face]>=0) && (bc[face]<world_size) &&
	(bc[face]!=world_rank);
      if( shared[face] ) range[face] = g->range[bc[face]]; 
    }

    // Begin receiving the particle counts

    for (face = 0; face < 6; face++) {
      if (shared[face]) {
	mp_size_recv_buffer( mp, f2b[face], sizeof(int) );
	mp_begin_recv( mp, f2b[face], sizeof(int), bc[face], f2rb[face] );
      }
    }
    
    // Load the particle send and local injection buffers
    do {

      particle_injector_t * RESTRICT ALIGNED(16) pi_send[6];

      // Presize the send and injection buffers
      //
      // Each buffer is large enough to hold one injector corresponding
      // to every mover in use (worst case, but plausible scenario in
      // beam simulations, is one buffer gets all the movers).
    
      int nm = 0;
      for (auto sp = vmprts.cbegin(); sp != vmprts.cend(); ++sp) {
	nm += sp->nm;
      }

      for( face=0; face<6; face++ )
	if( shared[face] ) {
	  mp_size_send_buffer( mp, f2b[face], 16+nm*sizeof(particle_injector_t) );
	  pi_send[face] = (particle_injector_t *)(((char *)mp_send_buffer(mp,f2b[face]))+16);
	  n_send[face] = 0;
	}

      if( max_ci<nm ) {
	particle_injector_t * new_ci = ci;
	FREE_ALIGNED( new_ci );
	MALLOC_ALIGNED( new_ci, nm, 16 );
	ci     = new_ci;
	max_ci = nm;
      }
      n_ci = 0;

      // For each species, load the movers

      for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
	const float   sp_q  = sp->q;
	const int32_t sp_id = sp->id;

	particle_t * RESTRICT ALIGNED(128) p0 = sp->p;
	int np = sp->np;

	particle_mover_t * RESTRICT ALIGNED(16)  pm = sp->pm + sp->nm - 1;
	nm = sp->nm;

	particle_injector_t * RESTRICT ALIGNED(16) pi;
	int i, voxel;
	int64_t nn;
      
	// Note that particle movers for each species are processed in
	// reverse order.  This allows us to backfill holes in the
	// particle list created by boundary conditions and/or
	// communication.  This assumes particle on the mover list are
	// monotonically increasing.  That is: pm[n].i > pm[n-1].i for
	// n=1...nm-1.  advance_p and inject_particle create movers with
	// property if all aged particle injection occurs after
	// advance_p and before this

	for( ; nm; pm--, nm-- ) {
	  i = pm->i;
	  voxel = p0[i].i;
	  face = voxel & 7;
	  voxel >>= 3;
	  p0[i].i = voxel;
	  nn = neighbor[ 6*voxel + face ];
        
	  // Absorb

	  if( nn==absorb_particles ) {
	    // Ideally, we would batch all rhob accumulations together
	    // for efficiency
	    accumulate_rhob(fa, p0+i, sp_q );
	    goto backfill;
	  }

	  // Send to a neighboring node

	  if( ((nn>=0) & (nn< rangel)) | ((nn>rangeh) & (nn<=rangem)) ) {
	    pi = &pi_send[face][n_send[face]++];
#         ifdef V4_ACCELERATION
	    copy_4x1( &pi->dx,    &p0[i].dx  );
	    copy_4x1( &pi->ux,    &p0[i].ux  );
	    copy_4x1( &pi->dispx, &pm->dispx );
#         else
	    pi->dx=p0[i].dx; pi->dy=p0[i].dy; pi->dz=p0[i].dz;
	    pi->i =nn - range[face];
	    pi->ux=p0[i].ux; pi->uy=p0[i].uy; pi->uz=p0[i].uz; pi->w=p0[i].w;
	    pi->dispx = pm->dispx; pi->dispy = pm->dispy; pi->dispz = pm->dispz;
	    pi->sp_id = sp_id;
#         endif 
	    (&pi->dx)[axis[face]] = dir[face];
	    pi->i                 = nn - range[face];
	    pi->sp_id             = sp_id;
	    goto backfill;
	  }

	  // Uh-oh: We fell through

	  WARNING(( "Unknown boundary interaction ... dropping particle "
		    "(species=%s)", sp->name ));

	backfill:

	  np--;
#       ifdef V4_ACCELERATION
	  copy_4x1( &p0[i].dx, &p0[np].dx );
	  copy_4x1( &p0[i].ux, &p0[np].ux );
#       else
	  p0[i] = p0[np];
#       endif

	}
      
	sp->np = np;
	sp->nm = 0;
      }

    } while(0);

    // Finish exchanging particle counts and start exchanging actual
    // particles.
  
    for (face=0; face<6; face++) {
      if (shared[face]) {
	*((int *)mp_send_buffer( mp, f2b[face] )) = n_send[face];
	mp_begin_send( mp, f2b[face], sizeof(int), bc[face], f2b[face] );
      }
    }
    for (face=0; face<6; face++) {
      if (shared[face])  {
	mp_end_recv( mp, f2b[face] );
	n_recv[face] = *((int *)mp_recv_buffer( mp, f2b[face] ));
	mp_size_recv_buffer( mp, f2b[face],
			     16+n_recv[face]*sizeof(particle_injector_t) );
	mp_begin_recv( mp, f2b[face], 16+n_recv[face]*sizeof(particle_injector_t),
		       bc[face], f2rb[face] );
      }
    }
    for (face=0; face<6; face++) {
      if (shared[face]) {
	mp_end_send( mp, f2b[face] );
	// FIXME: ASSUMES MP WON'T MUCK WITH REST OF SEND BUFFER. IF WE
	// DID MORE EFFICIENT MOVER ALLOCATION ABOVE, THIS WOULD BE
	// ROBUSTED AGAINST MP IMPLEMENTATION VAGARIES
	mp_begin_send( mp, f2b[face], 16+n_send[face]*sizeof(particle_injector_t),
		       bc[face], f2b[face] );
      }
    }

    do {
      // Unpack the species list for random acesss

      particle_t       * RESTRICT ALIGNED(32) sp_p[ MAX_SP];
      particle_mover_t * RESTRICT ALIGNED(32) sp_pm[MAX_SP];
      float sp_q[MAX_SP];
      int sp_np[MAX_SP];
      int sp_nm[MAX_SP];
      int sp_max_np[MAX_SP], n_dropped_particles[MAX_SP];
      int sp_max_nm[MAX_SP], n_dropped_movers[MAX_SP];

      if (vmprts.getNumSpecies() > MAX_SP) {
	ERROR(( "Update this to support more species" ));
      }
      for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
	sp_p[  sp->id ] = sp->p;
	sp_pm[ sp->id ] = sp->pm;
	sp_q[  sp->id ] = sp->q;
	sp_np[ sp->id ] = sp->np;
	sp_nm[ sp->id ] = sp->nm;
	sp_max_np[sp->id]=sp->max_np; n_dropped_particles[sp->id]=0;
	sp_max_nm[sp->id]=sp->max_nm; n_dropped_movers[sp->id]=0;
      }

      // Inject particles.  We do custom local injection first to
      // increase message overlap opportunities.

      face = 5;
      do {
	/**/  particle_t          * RESTRICT ALIGNED(32) p;
	/**/  particle_mover_t    * RESTRICT ALIGNED(16) pm;
	const particle_injector_t * RESTRICT ALIGNED(16) pi;
	int np, nm, n, id;
  
	face++;
	if (face == 7) {
	  face = 0;
	}

	if (face == 6) {
	  pi = ci;
	  n = n_ci;
	} else if (shared[face]) {
	  mp_end_recv( mp, f2b[face] );
	  pi = (const particle_injector_t *)
	    (((char *)mp_recv_buffer(mp,f2b[face]))+16);
	  n  = n_recv[face];
	} else {
	  continue;
	}
        
	// Reverse order injection is done to reduce thrashing of the
	// particle list (particles are removed reverse order so the
	// overall impact of removal + injection is to keep injected
	// particles in order).
	//
	// WARNING: THIS TRUSTS THAT THE INJECTORS (INCLUDING THOSE
	// RECEIVED FROM OTHER NODES) HAVE VALID PARTICLE IDS.
  
	pi += n-1;
	for( ; n; pi--, n-- ) {
	  id = pi->sp_id;
	  p  = sp_p[id];  np = sp_np[id];
	  pm = sp_pm[id]; nm = sp_nm[id];

	  if( np>=sp_max_np[id] ) { n_dropped_particles[id]++; continue; }
#       ifdef V4_ACCELERATION
	  copy_4x1(  &p[np].dx,    &pi->dx    );
	  copy_4x1(  &p[np].ux,    &pi->ux    );
#       else
	  p[np].dx=pi->dx; p[np].dy=pi->dy; p[np].dz=pi->dz; p[np].i=pi->i;
	  p[np].ux=pi->ux; p[np].uy=pi->uy; p[np].uz=pi->uz; p[np].w=pi->w;
#       endif
	  sp_np[id] = np+1;

	  if( nm>=sp_max_nm[id] ) { n_dropped_movers[id]++;    continue; }
#       ifdef V4_ACCELERATION
	  copy_4x1( &pm[nm].dispx, &pi->dispx );
	  pm[nm].i = np;
#       else
	  pm[nm].dispx=pi->dispx; pm[nm].dispy=pi->dispy; pm[nm].dispz=pi->dispz;
	  pm[nm].i=np;
#       endif
	  sp_nm[id] = nm + move_p(p, pm+nm, acc_block, g, sp_q[id]);
	}
      } while (face != 5);
  
      for (auto sp = vmprts.begin(); sp != vmprts.end(); ++sp) {
	if( n_dropped_particles[sp->id] )
	  WARNING(( "Dropped %i particles from species \"%s\".  Use a larger "
		    "local particle allocation in your simulation setup for "
		    "this species on this node.",
		    n_dropped_particles[sp->id], sp->name ));
	if( n_dropped_movers[sp->id] )
	  WARNING(( "%i particles were not completed moved to their final "
		    "location this timestep for species \"%s\".  Use a larger "
		    "local particle mover buffer in your simulation setup "
		    "for this species on this node.",
		    n_dropped_movers[sp->id], sp->name ));
	sp->np=sp_np[sp->id];
	sp->nm=sp_nm[sp->id];
      }

    } while(0);
  
    for (face = 0; face < 6; face++) {
      if (shared[face]) {
	mp_end_send(mp,f2b[face]);
      }
    }
  }
  
  
  void boundary_p(particle_bc_t *pbc_list, Particles& vmprts, FieldArray& fa,
		  Accumulator& accumulator)
  {
    boundary_p_(pbc_list, vmprts, fa, accumulator[0]);
  }

  // ----------------------------------------------------------------------
  // accumulate_rho_p

  void accumulate_rho_p(FieldArray& fa, typename Particles::const_iterator sp)
  {
    const particle_t * RESTRICT ALIGNED(128) p = sp->p;

    const float q_8V = sp->q*sp->g->r8V;
    const int np = sp->np;
    const int sy = sp->g->sy;
    const int sz = sp->g->sz;

    float w0, w1, w2, w3, w4, w5, w6, w7, dz;

    int n, v;

    // Load the grid data

    for( n=0; n<np; n++ ) {

      // Load the particle data

      w0 = p[n].dx;
      w1 = p[n].dy;
      dz = p[n].dz;
      v  = p[n].i;
      w7 = p[n].w*q_8V;

      // Compute the trilinear weights
      // Though the PPE should have hardware fma/fmaf support, it was
      // measured to be more efficient _not_ to use it here.  (Maybe the
      // compiler isn't actually generating the assembly for it.

#   define FMA( x,y,z) ((z)+(x)*(y))
#   define FNMS(x,y,z) ((z)-(x)*(y))
      w6=FNMS(w0,w7,w7);                    // q(1-dx)
      w7=FMA( w0,w7,w7);                    // q(1+dx)
      w4=FNMS(w1,w6,w6); w5=FNMS(w1,w7,w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
      w6=FMA( w1,w6,w6); w7=FMA( w1,w7,w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
      w0=FNMS(dz,w4,w4); w1=FNMS(dz,w5,w5); w2=FNMS(dz,w6,w6); w3=FNMS(dz,w7,w7);
      w4=FMA( dz,w4,w4); w5=FMA( dz,w5,w5); w6=FMA( dz,w6,w6); w7=FMA( dz,w7,w7);
#   undef FNMS
#   undef FMA

      // Reduce the particle charge to rhof

      fa[v      ].rhof += w0; fa[v      +1].rhof += w1;
      fa[v   +sy].rhof += w2; fa[v   +sy+1].rhof += w3;
      fa[v+sz   ].rhof += w4; fa[v+sz   +1].rhof += w5;
      fa[v+sz+sy].rhof += w6; fa[v+sz+sy+1].rhof += w7;
    }
  }

  void accumulate_rho_p(Particles& vmprts, FieldArray &vmflds)
  {
    for (auto sp = vmprts.cbegin(); sp != vmprts.cend(); ++sp) {
      TIC accumulate_rho_p(vmflds, sp); TOC(accumulate_rho_p, 1);
    }
  }

  // ----------------------------------------------------------------------
  // accumulate_rhob

  void accumulate_rhob(FieldArray& fa, const particle_t* p, float qsp)
  {
    float w0 = p->dx, w1 = p->dy, w2, w3, w4, w5, w6, w7, dz = p->dz;
    int v = p->i, x, y, z, sy = fa.g->sy, sz = fa.g->sz;
    const int nx = fa.g->nx, ny = fa.g->ny, nz = fa.g->nz;
    w7 = (qsp * fa.g->r8V) * p->w;

    // Compute the trilinear weights
    // See note in rhof for why FMA and FNMS are done this way.

# define FMA( x,y,z) ((z)+(x)*(y))
# define FNMS(x,y,z) ((z)-(x)*(y))
    w6=FNMS(w0,w7,w7);                    // q(1-dx)
    w7=FMA( w0,w7,w7);                    // q(1+dx)
    w4=FNMS(w1,w6,w6); w5=FNMS(w1,w7,w7); // q(1-dx)(1-dy), q(1+dx)(1-dy)
    w6=FMA( w1,w6,w6); w7=FMA( w1,w7,w7); // q(1-dx)(1+dy), q(1+dx)(1+dy)
    w0=FNMS(dz,w4,w4); w1=FNMS(dz,w5,w5); w2=FNMS(dz,w6,w6); w3=FNMS(dz,w7,w7);
    w4=FMA( dz,w4,w4); w5=FMA( dz,w5,w5); w6=FMA( dz,w6,w6); w7=FMA( dz,w7,w7);
# undef FNMS
# undef FMA

    // Adjust the weights for a corrected local accumulation of rhob.
    // See note in synchronize_rho why we must do this for rhob and not
    // for rhof.

    x  = v;    z = x/sz;
    if (z==1 ) w0 += w0, w1 += w1, w2 += w2, w3 += w3;
    if (z==nz) w4 += w4, w5 += w5, w6 += w6, w7 += w7;
    x -= sz*z; y = x/sy;
    if (y==1 ) w0 += w0, w1 += w1, w4 += w4, w5 += w5;
    if (y==ny) w2 += w2, w3 += w3, w6 += w6, w7 += w7;
    x -= sy*y;
    if (x==1 ) w0 += w0, w2 += w2, w4 += w4, w6 += w6;
    if (x==nx) w1 += w1, w3 += w3, w5 += w5, w7 += w7;

    // Reduce the particle charge to rhob

    fa[v      ].rhob += w0; fa[v      +1].rhob += w1;
    fa[v   +sy].rhob += w2; fa[v   +sy+1].rhob += w3;
    fa[v+sz   ].rhob += w4; fa[v+sz   +1].rhob += w5;
    fa[v+sz+sy].rhob += w6; fa[v+sz+sy+1].rhob += w7;
  }

  //private:
  vpic_simulation *simulation_;
};



template<class ParticlesBase, class FA, class IA, class AA, class HA>
struct PscParticles : ParticlesBase
{
  typedef ParticlesBase Base;
  typedef PscParticles<Base, FA, IA, AA, HA> Particles;
  typedef FA FieldArray;
  typedef IA Interpolator;
  typedef AA Accumulator;
  typedef HA HydroArray;

  using Base::Base;

  // ----------------------------------------------------------------------
  // accumulate_hydro_p
  
  // accumulate_hydro_p adds the hydrodynamic fields associated with the
  // supplied particle_list to the hydro array.  Trilinear interpolation
  // is used.  hydro is known at the nodes at the same time as particle
  // positions. No effort is made to fix up edges of the computational
  // domain.  All particles on the list must be inbounds.  Note, the
  // hydro jx,jy,jz are for diagnostic purposes only; they are not
  // accumulated with a charge conserving algorithm.

  static void accumulate_hydro_p(HydroArray& ha, typename Particles::const_iterator sp,
				 /*const*/ Interpolator& IP)
  {
    float c, qsp, mspc, qdt_2mc, qdt_4mc2, r8V;
    int np, stride_10, stride_21, stride_43;

    float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
    float w0, w1, w2, w3, w4, w5, w6, w7, t;
    int i, n;

    const particle_t *p = sp->p;

    c        = sp->g->cvac;
    qsp      = sp->q;
    mspc     = sp->m*c;
    qdt_2mc  = (qsp*sp->g->dt)/(2*mspc);
    qdt_4mc2 = qdt_2mc / (2*c);
    r8V      = sp->g->r8V;

    np        = sp->np;
    stride_10 = (VOXEL(1,0,0, sp->g->nx,sp->g->ny,sp->g->nz) -
		 VOXEL(0,0,0, sp->g->nx,sp->g->ny,sp->g->nz));
    stride_21 = (VOXEL(0,1,0, sp->g->nx,sp->g->ny,sp->g->nz) -
		 VOXEL(1,0,0, sp->g->nx,sp->g->ny,sp->g->nz));
    stride_43 = (VOXEL(0,0,1, sp->g->nx,sp->g->ny,sp->g->nz) -
		 VOXEL(1,1,0, sp->g->nx,sp->g->ny,sp->g->nz));

    for( n=0; n<np; n++ ) {

      // Load the particle
      dx = p[n].dx;
      dy = p[n].dy;
      dz = p[n].dz;
      i  = p[n].i;
      ux = p[n].ux;
      uy = p[n].uy;
      uz = p[n].uz;
      w  = p[n].w;
    
      // Half advance E
      ux += qdt_2mc*((IP[i].ex+dy*IP[i].dexdy) + dz*(IP[i].dexdz+dy*IP[i].d2exdydz));
      uy += qdt_2mc*((IP[i].ey+dz*IP[i].deydz) + dx*(IP[i].deydx+dz*IP[i].d2eydzdx));
      uz += qdt_2mc*((IP[i].ez+dx*IP[i].dezdx) + dy*(IP[i].dezdy+dx*IP[i].d2ezdxdy));

      // Boris rotation - Interpolate B field
      w5 = IP[i].cbx + dx*IP[i].dcbxdx;
      w6 = IP[i].cby + dy*IP[i].dcbydy;
      w7 = IP[i].cbz + dz*IP[i].dcbzdz;

      // Boris rotation - curl scalars (0.5 in v0 for half rotate) and
      // kinetic energy computation. Note: gamma-1 = |u|^2 / (gamma+1)
      // is the numerically accurate way to compute gamma-1
      ke_mc = ux*ux + uy*uy + uz*uz; // ke_mc = |u|^2 (invariant)
      vz = sqrt(1+ke_mc);            // vz = gamma    (invariant)
      ke_mc *= c/(vz+1);             // ke_mc = c|u|^2/(gamma+1) = c*(gamma-1)
      vz = c/vz;                     // vz = c/gamma
      w0 = qdt_4mc2*vz;
      w1 = w5*w5 + w6*w6 + w7*w7;    // |cB|^2
      w2 = w0*w0*w1;
      w3 = w0*(1+(1./3.)*w2*(1+0.4*w2));
      w4 = w3/(1 + w1*w3*w3); w4 += w4;

      // Boris rotation - uprime
      w0 = ux + w3*( uy*w7 - uz*w6 );
      w1 = uy + w3*( uz*w5 - ux*w7 );
      w2 = uz + w3*( ux*w6 - uy*w5 );

      // Boris rotation - u
      ux += w4*( w1*w7 - w2*w6 );
      uy += w4*( w2*w5 - w0*w7 );
      uz += w4*( w0*w6 - w1*w5 );

      // Compute physical velocities
      vx  = ux*vz;
      vy  = uy*vz;
      vz *= uz;

      // Compute the trilinear coefficients
      w0  = r8V*w;    // w0 = (1/8)(w/V)
      dx *= w0;       // dx = (1/8)(w/V) x
      w1  = w0+dx;    // w1 = (1/8)(w/V) + (1/8)(w/V)x = (1/8)(w/V)(1+x)
      w0 -= dx;       // w0 = (1/8)(w/V) - (1/8)(w/V)x = (1/8)(w/V)(1-x)
      w3  = 1+dy;     // w3 = 1+y
      w2  = w0*w3;    // w2 = (1/8)(w/V)(1-x)(1+y)
      w3 *= w1;       // w3 = (1/8)(w/V)(1+x)(1+y)
      dy  = 1-dy;     // dy = 1-y
      w0 *= dy;       // w0 = (1/8)(w/V)(1-x)(1-y)
      w1 *= dy;       // w1 = (1/8)(w/V)(1+x)(1-y)
      w7  = 1+dz;     // w7 = 1+z
      w4  = w0*w7;    // w4 = (1/8)(w/V)(1-x)(1-y)(1+z) = (w/V) trilin_0 *Done
      w5  = w1*w7;    // w5 = (1/8)(w/V)(1+x)(1-y)(1+z) = (w/V) trilin_1 *Done
      w6  = w2*w7;    // w6 = (1/8)(w/V)(1-x)(1+y)(1+z) = (w/V) trilin_2 *Done
      w7 *= w3;       // w7 = (1/8)(w/V)(1+x)(1+y)(1+z) = (w/V) trilin_3 *Done
      dz  = 1-dz;     // dz = 1-z
      w0 *= dz;       // w0 = (1/8)(w/V)(1-x)(1-y)(1-z) = (w/V) trilin_4 *Done
      w1 *= dz;       // w1 = (1/8)(w/V)(1+x)(1-y)(1-z) = (w/V) trilin_5 *Done
      w2 *= dz;       // w2 = (1/8)(w/V)(1-x)(1+y)(1-z) = (w/V) trilin_6 *Done
      w3 *= dz;       // w3 = (1/8)(w/V)(1+x)(1+y)(1-z) = (w/V) trilin_7 *Done

      // Accumulate the hydro fields
#   define ACCUM_HYDRO( wn)						\
      t  = qsp*wn;        /* t  = (qsp w/V) trilin_n */			\
      ha[i].jx  += t*vx;						\
      ha[i].jy  += t*vy;						\
      ha[i].jz  += t*vz;						\
      ha[i].rho += t;							\
      t  = mspc*wn;       /* t = (msp c w/V) trilin_n */		\
      dx = t*ux;          /* dx = (px w/V) trilin_n */			\
      dy = t*uy;							\
      dz = t*uz;							\
      ha[i].px  += dx;							\
      ha[i].py  += dy;							\
      ha[i].pz  += dz;							\
      ha[i].ke  += t*ke_mc;						\
      ha[i].txx += dx*vx;						\
      ha[i].tyy += dy*vy;						\
      ha[i].tzz += dz*vz;						\
      ha[i].tyz += dy*vz;						\
      ha[i].tzx += dz*vx;						\
      ha[i].txy += dx*vy

      /**/            ACCUM_HYDRO(w0); // Cell i,j,k
      i += stride_10; ACCUM_HYDRO(w1); // Cell i+1,j,k
      i += stride_21; ACCUM_HYDRO(w2); // Cell i,j+1,k
      i += stride_10; ACCUM_HYDRO(w3); // Cell i+1,j+1,k
      i += stride_43; ACCUM_HYDRO(w4); // Cell i,j,k+1
      i += stride_10; ACCUM_HYDRO(w5); // Cell i+1,j,k+1
      i += stride_21; ACCUM_HYDRO(w6); // Cell i,j+1,k+1
      i += stride_10; ACCUM_HYDRO(w7); // Cell i+1,j+1,k+1

#   undef ACCUM_HYDRO
    }
  }

  // ----------------------------------------------------------------------
  // uncenter_p
  
  static void uncenter_p_pipeline(species_t *sp, /*const*/ Interpolator& interpolator,
				  int off, int cnt)
  {
    const interpolator_t*  f;
    // For backward half advance
    const float qdt_2mc = -(sp->q * sp->g->dt) / (2*sp->m * sp->g->cvac);
    const float qdt_4mc = 0.5 * qdt_2mc; // For backward half rotate
    const float one = 1.;
    const float one_third = 1./3.;
    const float two_fifteenths = 2./15.;

    float dx, dy, dz, ux, uy, uz;
    float hax, hay, haz, cbx, cby, cbz;
    float v0, v1, v2, v3, v4;

    int ii;

    int n = cnt;
    particle_t *p = sp->p + off;

    // Process particles for this pipeline

    for (; n; n--, p++) {
      dx   = p->dx;                            // Load position
      dy   = p->dy;
      dz   = p->dz;
      ii   = p->i;
      f    = &interpolator[ii];                // Interpolate E
      hax  = qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
		       dz*( f->dexdz + dy*f->d2exdydz ) );
      hay  = qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
		       dx*( f->deydx + dz*f->d2eydzdx ) );
      haz  = qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
		       dy*( f->dezdy + dx*f->d2ezdxdy ) );
      cbx  = f->cbx + dx*f->dcbxdx;            // Interpolate B
      cby  = f->cby + dy*f->dcbydy;
      cbz  = f->cbz + dz*f->dcbzdz;
      ux   = p->ux;                            // Load momentum
      uy   = p->uy;
      uz   = p->uz;
      v0   = qdt_4mc/(float)sqrt(one + (ux*ux + (uy*uy + uz*uz)));
      /**/                                     // Boris - scalars
      v1   = cbx*cbx + (cby*cby + cbz*cbz);
      v2   = (v0*v0)*v1;
      v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
      v4   = v3/(one+v1*(v3*v3));
      v4  += v4;
      v0   = ux + v3*( uy*cbz - uz*cby );      // Boris - uprime
      v1   = uy + v3*( uz*cbx - ux*cbz );
      v2   = uz + v3*( ux*cby - uy*cbx );
      ux  += v4*( v1*cbz - v2*cby );           // Boris - rotation
      uy  += v4*( v2*cbx - v0*cbz );
      uz  += v4*( v0*cby - v1*cbx );
      ux  += hax;                              // Half advance E
      uy  += hay;
      uz  += haz;
      p->ux = ux;                              // Store momentum
      p->uy = uy;
      p->uz = uz;
    }
  }

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

  static void uncenter_p_pipeline_v4(species_t *sp, /*const*/ Interpolator& interpolator,
				     int off, int cnt)
  {
    using namespace v4;
    const interpolator_t * ALIGNED(128) f0 = interpolator.data();

    particle_t           * ALIGNED(128) p;
    const float          * ALIGNED(16)  vp0;
    const float          * ALIGNED(16)  vp1;
    const float          * ALIGNED(16)  vp2;
    const float          * ALIGNED(16)  vp3;

    const float _qdt_2mc = (sp->q * sp->g->dt) / (2*sp->m * sp->g->cvac);
    
    const v4float qdt_2mc(    -_qdt_2mc); // For backward half advance
    const v4float qdt_4mc(-0.5*_qdt_2mc); // For backward half Boris rotate
    const v4float one(1.);
    const v4float one_third(1./3.);
    const v4float two_fifteenths(2./15.);

    v4float dx, dy, dz, ux, uy, uz, q;
    v4float hax, hay, haz, cbx, cby, cbz;
    v4float v0, v1, v2, v3, v4, v5;
    v4int ii;

    // Determine which particle quads this pipeline processes

    assert(off % 4 == 0);
    assert(cnt % 4 == 0);
    int nq = cnt >> 2;
    p = sp->p + off;

    // Process the particle quads for this pipeline

    for( ; nq; nq--, p+=4 ) {
      load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

      // Interpolate fields
      vp0 = (const float * ALIGNED(16))(f0 + ii(0));
      vp1 = (const float * ALIGNED(16))(f0 + ii(1));
      vp2 = (const float * ALIGNED(16))(f0 + ii(2));
      vp3 = (const float * ALIGNED(16))(f0 + ii(3));
      load_4x4_tr(vp0,  vp1,  vp2,  vp3,  hax,v0,v1,v2); hax = qdt_2mc*fma( fma( dy, v2, v1 ), dz, fma( dy, v0, hax ) );
      load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,hay,v3,v4,v5); hay = qdt_2mc*fma( fma( dz, v5, v4 ), dx, fma( dz, v3, hay ) );
      load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,haz,v0,v1,v2); haz = qdt_2mc*fma( fma( dx, v2, v1 ), dy, fma( dx, v0, haz ) );
      load_4x4_tr(vp0+12,vp1+12,vp2+12,vp3+12,cbx,v3,cby,v4); cbx = fma( v3, dx, cbx );
      /**/                                                    cby = fma( v4, dy, cby );
      load_4x2_tr(vp0+16,vp1+16,vp2+16,vp3+16,cbz,v5);        cbz = fma( v5, dz, cbz );

      // Update momentum
      load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,ux,uy,uz,q);
      /**/                                              // Could use load_4x3_tr
      v0  = qdt_4mc*rsqrt( one + fma( ux,ux, fma( uy,uy, uz*uz ) ) );
      v1  = fma( cbx,cbx, fma( cby,cby, cbz*cbz ) );
      v2  = (v0*v0)*v1;
      v3  = v0*fma( v2, fma( v2, two_fifteenths, one_third ), one );
      v4  = v3*rcp( fma( v3*v3, v1, one ) ); v4 += v4;
      v0  = fma( fms( uy,cbz, uz*cby ), v3, ux );
      v1  = fma( fms( uz,cbx, ux*cbz ), v3, uy );
      v2  = fma( fms( ux,cby, uy*cbx ), v3, uz );
      ux  = fma( fms( v1,cbz, v2*cby ), v4, ux );
      uy  = fma( fms( v2,cbx, v0*cbz ), v4, uy );
      uz  = fma( fms( v0,cby, v1*cbx ), v4, uz );
      ux += hax;
      uy += hay;
      uz += haz;
      store_4x4_tr(ux,uy,uz,q,&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux);
      /**/                                              // Could use store_4x3_tr
    }
  }

#endif
  
  static void uncenter_p(species_t *sp, Interpolator& interpolator)
  {
    //    ::uncenter_p(sp, &interpolator); return;
    int cnt = sp->np & ~15;
#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)
    uncenter_p_pipeline_v4(sp, interpolator, 0, cnt);
#else
    uncenter_p_pipeline(sp, interpolator, 0, cnt);
#endif
    uncenter_p_pipeline(sp, interpolator, cnt, sp->np - cnt);
  }

  // ----------------------------------------------------------------------
  // sort_p

  static void sort_p(species_t * sp)
  {
    sp->last_sorted = sp->g->step;

    int n_prts = sp->np;
    int vl = VOXEL(1,1,1,                         sp->g->nx,sp->g->ny,sp->g->nz);
    int vh = VOXEL(sp->g->nx,sp->g->ny,sp->g->nz, sp->g->nx,sp->g->ny,sp->g->nz) + 1;

    static int * RESTRICT ALIGNED(128) next;
    if (!next) {
      MALLOC_ALIGNED(next, sp->g->nv * sizeof(*next), 128);
    }
    int * RESTRICT ALIGNED(128) partition = sp->partition;

    static particle_t * RESTRICT ALIGNED(128) p_aux;
    static size_t n_alloced;
    if (n_prts > n_alloced) {
      FREE_ALIGNED(p_aux);
      MALLOC_ALIGNED(p_aux, n_prts * sizeof(*p_aux), 128);
      n_alloced = n_prts;
    }
    particle_t * RESTRICT ALIGNED(128) p = sp->p;

    // zero counts
    CLEAR(&next[vl], vh - vl);
    
    // find counts
    for (int i = 0; i < n_prts; i++) {
      next[p[i].i]++;
    }
    
    // prefix sum
    int sum = 0;
    for (int v = vl; v < vh; v++) {
      int count = next[v];
      next[v] = sum;
      partition[v] = sum;
      sum += count;
    }
    partition[vh] = sum;
    
    // reorder
    for(int i = 0; i < n_prts; i++) {
      int v = p[i].i;
      int j = next[v]++;
      p_aux[j] = p[i];
    }

    // fix up unused part of partition
    CLEAR(partition, vl);
    for (int i = vh; i < sp->g->nv; i++) {
      partition[i] = n_prts;
    }

    // OPT: just swap pointer?
    COPY(p, p_aux, n_prts);
  }

};



#endif

