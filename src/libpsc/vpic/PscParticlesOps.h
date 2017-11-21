
#ifndef PSC_PARTICLES_OPS
#define PSC_PARTICLES_OPS

template<class P>
struct PscParticlesOps {
  typedef P Particles;
  
  PscParticlesOps(vpic_simulation *simulation) : simulation_(simulation) { }

  void inject_particle(Particles *vmprts, int patch, const struct psc_particle_inject *prt)
  {
    assert(patch == 0);
    assert(simulation_->accumulator_array);
    species_t *sp = &*vmprts->find_id(prt->kind); //::find_species_id(prt->kind, vmprts->sl_);
    assert(sp);

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

    if (update_rhob) accumulate_rhob( simulation_->field_array->f, p, grid, -sp->q );

    if (age!=0) {
      if( sp->nm>=sp->max_nm )
	WARNING(( "No movers available to age injected  particle" ));
      particle_mover_t * pm = sp->pm + sp->nm;
      age *= grid->cvac*grid->dt/sqrt( ux*ux + uy*uy + uz*uz + 1 );
      pm->dispx = ux*age*grid->rdx;
      pm->dispy = uy*age*grid->rdy;
      pm->dispz = uz*age*grid->rdz;
      pm->i     = sp->np-1;
      sp->nm += move_p( sp->p, pm, simulation_->accumulator_array->a, grid, sp->q );
    }
    
  }

private:
  vpic_simulation *simulation_;
};

#endif

