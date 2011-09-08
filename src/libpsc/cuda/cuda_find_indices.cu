
#include <psc_cuda.h>

// FIXME, specific to 1x8x8, should be in ! .cu, so that cell_map works

EXTERN_C void
create_indices_host(unsigned int *cnis, struct cell_map *map,
		    particles_cuda_t *pp, struct psc_patch *patch)
{
  for (int i = 0; i < pp->n_part; i++) {
    particles_cuda_dev_t *p = &pp->h_part;
    particle_cuda_real_t dxi = 1.f / ppsc->dx[0];
    particle_cuda_real_t dyi = 1.f / ppsc->dx[1];
    particle_cuda_real_t dzi = 1.f / ppsc->dx[2];
    particle_cuda_real_t xi[3] = {
      (p->xi4[i].x - patch->xb[0]) * dxi,
      (p->xi4[i].y - patch->xb[1]) * dyi,
      (p->xi4[i].z - patch->xb[2]) * dzi };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = particle_base_real_nint(xi[d]);
    }
    
    int idx = (((pos[2] / 8) * (patch->ldims[1] / 8) + (pos[1] / 8)) << 6) |
      ((pos[2] & 4) << 3) |
      ((pos[2] & 2) << 2) |
      ((pos[2] & 1) << 1) |
      ((pos[1] & 4) << 2) |
      ((pos[1] & 2) << 1) |
      ((pos[1] & 1) << 0);
    cnis[i] = idx;

    assert(cnis[i] < map->N);
  }
}

EXTERN_C void
particles_cuda_copy_to_device(particles_cuda_t *pp)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *h_part = &pp->h_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  check(cudaMemcpy(d_part->xi4, h_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_part->pxi4, h_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
}

