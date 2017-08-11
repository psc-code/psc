
static void
exchange_particles_pre(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts, int p)
{
  struct psc *psc = bnd->psc;
  struct ddc_particles *ddcp = bnd->ddcp;

  int n_send = get_n_send(mprts, p);
  struct psc_patch *patch = &psc->patch[p];
  particle_real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = patch->ldims[d] * patch->dx[d];
  }
  particle_real_t *b_dxi = get_b_dxi(mprts, p);
  int *b_mx = get_b_mx(mprts, p);
  
  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  struct ddcp_patch *ddcp_patch = &ddcp->patches[p];
  ddcp_patch->head = get_head(mprts, p);
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    particle_buf_resize(&ddcp_patch->nei[dir1].send_buf, 0);
  }
  int n_end = ddcp_patch->head + n_send;
  for (int n = ddcp_patch->head; n < n_end; n++) {
    particle_t *prt = xchg_get_one(mprts, p, n);
    particle_real_t *xi = &prt->xi;
    particle_real_t *pxi = &prt->pxi;
    
    bool drop = false;
    int dir[3];
    for (int d = 0; d < 3; d++) {
      int bi = particle_real_fint(xi[d] * b_dxi[d]);
      if (bi < 0) {
	// FIXME, assumes every patch has same dimensions
	if (patch->off[d] != 0 || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	  xi[d] += xm[d];
	  dir[d] = -1;
	  bi = particle_real_fint(xi[d] * b_dxi[d]);
	  if (bi >= b_mx[d]) {
	    xi[d] = 0.;
	    dir[d] = 0;
	  }
	} else {
	  switch (psc->domain.bnd_part_lo[d]) {
	  case BND_PART_REFLECTING:
	    xi[d]  = -xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    break;
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else if (bi >= b_mx[d]) {
	if (patch->off[d] + patch->ldims[d] != psc->domain.gdims[d] ||
	    psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	  xi[d] -= xm[d];
	  dir[d] = +1;
	  bi = particle_real_fint(xi[d] * b_dxi[d]);
	  if (bi < 0) {
	    xi[d] = 0.;
	  }
	} else {
	  switch (psc->domain.bnd_part_hi[d]) {
	  case BND_PART_REFLECTING:
	    xi[d] = 2.f * xm[d] - xi[d];
	    pxi[d] = -pxi[d];
	    dir[d] = 0;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] *= (1. - 1e-6);
	    }
	    break;
	  case BND_PART_ABSORBING:
	    drop = true;
	    break;
	  default:
	    assert(0);
	  }
	}
      } else {
	dir[d] = 0;
      }
      if (!drop) {
	if (xi[d] < 0.f && xi[d] > -1e-6f) {
	  mprintf("d %d xi %g\n", d, xi[d]);
	  xi[d] = 0.f;
	}
	assert(xi[d] >= 0.f);
	assert(xi[d] <= xm[d]);
      }
    }
    if (!drop) {
      if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	xchg_append(mprts, p, ddcp_patch, prt);
      } else {
	ddc_particles_queue(ddcp, ddcp_patch, dir, prt);
      }
    }
  }
}

