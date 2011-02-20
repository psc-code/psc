
#include "output_fields.h"

// ----------------------------------------------------------------------
// copy_to_global helper

static void
copy_to_global(fields_base_real_t *fld, fields_base_real_t *buf,
	       int *ilo, int *ihi, int *ilg, int *img)
{
  int *gdims = psc.domain.gdims;
  int my = gdims[1];
  int mx = gdims[0];

  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	fld[(iz * my + iy) * mx + ix] =
	  buf[((iz - ilg[2]) * img[1] + iy - ilg[1]) * img[0] + ix - ilg[0]];
      }
    }
  }
}

// ----------------------------------------------------------------------
// write_fields_combine

void
write_fields_combine(struct psc_fields_list *list, 
		     void (*write_field)(void *ctx, fields_base_t *fld), void *ctx)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  foreach_patch(patch) {
    for (int m = 0; m < list->nr_flds; m++) {
      int s_ilo[3], s_ihi[3], s_ilg[3], s_img[3];
      fields_base_real_t *s_data = &F3_BASE(&list->flds[m], 0,
					    -psc.ibn[0], -psc.ibn[1], -psc.ibn[2]);
      
      for (int d = 0; d < 3; d++) {
	s_ilo[d] = psc.patch[patch].off[d];
	s_ihi[d] = psc.patch[patch].off[d] + psc.patch[patch].ldims[d];
	s_ilg[d] = psc.patch[patch].off[d] - psc.ibn[d];
	s_img[d] = psc.patch[patch].ldims[d] + 2 * psc.ibn[d];
      }
      
      if (rank != 0) {
	MPI_Send(s_ilo, 3, MPI_INT, 0, 100, MPI_COMM_WORLD);
	MPI_Send(s_ihi, 3, MPI_INT, 0, 101, MPI_COMM_WORLD);
	MPI_Send(s_ilg, 3, MPI_INT, 0, 102, MPI_COMM_WORLD);
	MPI_Send(s_img, 3, MPI_INT, 0, 103, MPI_COMM_WORLD);
	unsigned int sz = fields_base_size(&list->flds[m]);
	MPI_Send(s_data, sz, MPI_FIELDS_BASE_REAL, 0, 104, MPI_COMM_WORLD);
      } else { // rank == 0
	fields_base_t fld;
	fields_base_alloc(&fld, (int []) { 0, 0, 0}, psc.domain.gdims, 1);
	fld.name[0] = strdup(list->flds[m].name[0]);
	
	for (int n = 0; n < size; n++) {
	  int ilo[3], ihi[3], ilg[3], img[3];
	  fields_base_real_t *buf;
	  
	  if (n == 0) {
	    for (int d = 0; d < 3; d++) {
	      ilo[d] = s_ilo[d];
	      ihi[d] = s_ihi[d];
	      ilg[d] = s_ilg[d];
	      img[d] = s_img[d];
	    }
	    buf = s_data;
	  } else {
	    MPI_Recv(ilo, 3, MPI_INT, n, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(ihi, 3, MPI_INT, n, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(ilg, 3, MPI_INT, n, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(img, 3, MPI_INT, n, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    int ntot = img[0] * img[1] * img[2];
	    buf = calloc(ntot, sizeof(*buf));
	    MPI_Recv(buf, ntot, MPI_FIELDS_BASE_REAL, n, 104, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	  }
	  /* printf("[%d] ilo %d %d %d ihi %d %d %d\n", rank, ilo[0], ilo[1], ilo[2], */
	  /*        ihi[0], ihi[1], ihi[2]); */
	  copy_to_global(&F3_BASE(&fld, 0, 0,0,0), buf, ilo, ihi, ilg, img);
	  if (n != 0) {
	    free(buf);
	  }
	}
	write_field(ctx, &fld);
	fields_base_free(&fld);
      }
    }
  }
}

