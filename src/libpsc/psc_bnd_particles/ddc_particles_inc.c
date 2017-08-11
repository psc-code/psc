
#include <mrc_ddc.h>

#include <string.h>

static void
ddc_particles_queue(struct ddc_particles *ddcp, struct ddcp_patch *patch,
		    int dir[3], void *p)
{
  struct ddcp_nei *nei = &patch->nei[mrc_ddc_dir2idx(dir)];

  if (nei->n_send == nei->send_buf_size) {
    // reallocate a larger buffer, doubling buffer size each time
    assert(nei->send_buf_size > 0);
    nei->send_buf_size *= 2;
    nei->send_buf = realloc(nei->send_buf, 
			    nei->send_buf_size * ddcp->size_of_particle);
  }
  memcpy(nei->send_buf + nei->n_send * ddcp->size_of_particle, p,
	 ddcp->size_of_particle);
  nei->n_send++;
}

