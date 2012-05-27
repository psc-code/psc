/*
 *  XDMF/HDF5 particle output by Simon Jagoda
 *  Compact version, 1 pair of files per node (total, not per step)
 */



#include "psc_output_particles_private.h"
#include "psc_output_particles_xdmf.h"
#include "mrc_io.h"
#include "mpi.h"
#include "psc.h"
#include "mrc_domain.h"
#include "mrc_fld.h"
#include <mrc_params.h>
#include <assert.h>
#include <string.h>
#include "hdf5.h"

#include <mrc_profile.h>

#define to_psc_output_particles_xdmf_compact(out) ((struct psc_output_particles_xdmf_compact *)((out)->obj.subctx))
#define to_psc_output_particles_xdmf(out) ((struct psc_output_particles_xdmf *)((out)->obj.subctx))

//----------------------------------------------------------------------
// static function declarations now here
// because the compiler dislikes having them in the header
static void create_hdf5 (struct psc_output_particles_xdmf *out_c);
static void create_xdmf (struct psc_output_particles_xdmf *out_c);
static void write_node_hdf5(struct psc_output_particles_xdmf *out_c, mparticles_base_t *particles_in);
static void write_node_xdmf(struct psc_output_particles_xdmf *out_c, mparticles_base_t *particles_in);
static void write_patch_xdmf(struct psc_output_particles_xdmf *out_c, int *partnumber, int *patch);
static void write_patch_hdf5(struct psc_output_particles_xdmf *out_c, struct psc_particles *patchdata, int *patch);
static void close_xdmf(struct psc_output_particles_xdmf *out_c);
static void close_hdf5(struct psc_output_particles_xdmf *out_c); 

//----------------------------------------------------------------------
// xdmf master file linking the node xdmfs
static void create_master_xdmf(struct psc_output_particles_xdmf *out_c) {
	char filename[30+strlen(out_c->data_dir)];
	sprintf(filename,"%s/particles_master.xdmf",out_c->data_dir);
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    FILE *master = fopen (filename,"w");
	fprintf(master,"<?xml version=\"1.0\" ?> \n");
	fprintf(master,"<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'> \n");
	fprintf(master,"<Domain> \n"); 
	fprintf(master,"<Grid GridType=\"Collection\" CollectionType=\"Spatial\"> \n");
	for(int node=0; node<mpi_size; node++) {
		fprintf(master,"<xi:include href=\"./%s.n%06d.xdmf\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\"/> \n",out_c->basename,node);
	}
	fprintf(master,"</Grid> \n");
	fprintf(master,"</Domain> \n");
	fprintf(master,"</Xdmf> \n");
	fclose(master);
}

//----------------------------------------------------------------------
// set some variables -> call once at t=0
static void initialize_output (struct psc_output_particles_xdmf *out_c) {
	out_c->typenames[0] = "x";
    out_c->typenames[1] = "y";
    out_c->typenames[2] = "z"; 
    out_c->typenames[3] = "px";
    out_c->typenames[4] = "py";	//names hardcoded for now
    out_c->typenames[5] = "pz";
    out_c->typenames[6] = "q";
    out_c->typenames[7] = "m";
    out_c->typenames[8] = "w";
    out_c->basename = "particles";
    out_c->write_particles=(bool*) malloc(ppsc->prm.nmax*sizeof(bool));

}


//----------------------------------------------------------------------
// create files for output
void create_hdf5 (struct psc_output_particles_xdmf *out_c) {
    char charbuffer[30+strlen(out_c->data_dir)+strlen(out_c->basename)]; 
    int node;
    MPI_Comm_rank(MPI_COMM_WORLD, &node);	
    sprintf(out_c->hdf5filename,"%s.n%06d.h5",out_c->basename,node);
    sprintf(charbuffer,"%s/%s",out_c->data_dir,out_c->hdf5filename);
    out_c->hdf5 = H5Fcreate (charbuffer, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
}
static void create_xdmf(struct psc_output_particles_xdmf *out_c) {
    char filename[30+strlen(out_c->data_dir)+strlen(out_c->basename)];  
    int node;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &node);
    sprintf(filename,"%s/%s.n%06d.xdmf",out_c->data_dir,out_c->basename,node);
    out_c->xdmf = fopen (filename,"w");
	fprintf(out_c->xdmf,"<?xml version=\"1.0\" ?> \n");
	fprintf(out_c->xdmf,"<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'> \n");
	fprintf(out_c->xdmf,"<Domain> \n");  
	fprintf(out_c->xdmf,"<Grid GridType=\"Collection\" CollectionType=\"Temporal\"> \n"); 
	fprintf(out_c->xdmf,"<Time TimeType=\"List\"> \n");
	fprintf(out_c->xdmf,"<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"%i\"> \n",out_c->times_written);
	for(int i=0; i<=ppsc->timestep; i++) {
		if(out_c->write_particles[i]) {
			fprintf(out_c->xdmf,"%e \n",ppsc->dt*i);
		}
	} 
	fprintf(out_c->xdmf,"</DataItem> \n");
	fprintf(out_c->xdmf,"</Time> \n");
	
}

//---------------------------------------------------------------------
//write container and data for the whole node
static void write_node_xdmf(struct psc_output_particles_xdmf *out_c, mparticles_base_t *particles_in) {
	fprintf(out_c->xdmf,"<Grid GridType=\"Collection\" CollectionType=\"Spatial\"> \n");
	for(int i=0; i<particles_in->nr_patches; i++) {
	  struct psc_particles *prts_base = psc_mparticles_get_patch(particles_in, i);
	  if(prts_base->n_part>0) {write_patch_xdmf(out_c, &(prts_base->n_part), &i);}
	}
	fprintf(out_c->xdmf,"</Grid>"); 
}
static void write_node_hdf5(struct psc_output_particles_xdmf *out_c, mparticles_base_t *particles_in)
{  
  for (int p = 0; p < particles_in->nr_patches; p++) {	
    struct psc_particles *prts_base = psc_mparticles_get_patch(particles_in, p);
    struct psc_particles *prts = psc_particles_get_as(prts_base, "c", p);
    if(prts->n_part>0) {write_patch_hdf5(out_c, prts, &p);}    
    psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
  } 
}


//----------------------------------------------------------------------
//write container and data for one patch
static void write_patch_xdmf(struct psc_output_particles_xdmf *out_c, int *partnumber, int *patch) {
	fprintf(out_c->xdmf,"<Grid GridType=\"Uniform\"> \n ");
	fprintf(out_c->xdmf,"<Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%i 1 1\" /> \n",*partnumber);
	fprintf(out_c->xdmf,"<Geometry GeometryType=\"Origin_DxDyDz\"> \n");
	fprintf(out_c->xdmf,"<DataItem Name=\"Origin\" Datatype=\"Float\" Dimensions=\"3\" Format=\"XML\"> \n");
	fprintf(out_c->xdmf,"0 0 0 \n");
	fprintf(out_c->xdmf,"</DataItem> \n");
	fprintf(out_c->xdmf,"<DataItem Name=\"DxDyDz\" Datatype=\"Float\" Dimensions=\"3\" Format=\"XML\"> \n");
	fprintf(out_c->xdmf,"1 1 1 \n");
	fprintf(out_c->xdmf,"</DataItem> \n");
	fprintf(out_c->xdmf,"</Geometry> \n");
	for(int offset=0; offset<9; offset++) {
		if(!out_c->write_variable[offset]) {continue;}	//skip this iteration if variable should not be written  
		fprintf(out_c->xdmf,"<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\"> \n",out_c->typenames[offset]);
		fprintf(out_c->xdmf,"<DataItem ItemType=\"Uniform\" Dimensions=\"%i 1 1\" NumberType=\"Float\" Format=\"HDF\"> \n", *partnumber);
		fprintf(out_c->xdmf,"./%s:/%s.%06d.%06d \n",out_c->hdf5filename,out_c->typenames[offset],ppsc->timestep,*patch);
		fprintf(out_c->xdmf,"</DataItem> \n");
		fprintf(out_c->xdmf,"</Attribute> \n");
	}
	fprintf(out_c->xdmf,"</Grid> \n"); 
}
static void write_patch_hdf5(struct psc_output_particles_xdmf *out_c, struct psc_particles *patchdata, int *patch) {    
    hid_t       memspace, filespace, dataset, dcpl;                                         
    hsize_t     memdims[1] = {9*patchdata->n_part}; 
    hsize_t	filedims[1] = {patchdata->n_part};
    hsize_t	start[1],
                stride[1],
                count[1],
                block[1];                     
    char datasetname[80];
  
    struct psc_particles_c *c = psc_particles_c(patchdata);
    //create dataspace for memory and file
    //these wont change for different particle properties
    memspace = H5Screate_simple (1, memdims, NULL);
    filespace = H5Screate_simple (1, filedims, NULL);
   
    //dataset creation property list (tm)
    //chunksize = everything there is
    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk (dcpl, 1, filedims);

    //constant hyperslab properties	
    stride[0] = 9;	//take each n-th element
    count[0] = patchdata->n_part;	//number of repetitions: array has to be of size stride*count or larger
    block[0] = 1;	//best not to touch this
                
    //loop over all particle variables                           
    for(int offset=0; offset<9; offset++) {    
    	if(!out_c->write_variable[offset]) {continue;}	//skip this iteration if variable should not be written              
            
     	sprintf(datasetname,"%s.%06d.%06d",out_c->typenames[offset],ppsc->timestep,*patch);
    	//hyperslab offset: determines particle variable to be written
    	start[0] = offset;	
    
    	//set hyperslab
    	H5Sselect_hyperslab (memspace, H5S_SELECT_SET, start, stride, count, block);

    	//create dataset: have to do this once for every variable, since name changes
    	dataset = H5Dcreate2 (out_c->hdf5, datasetname, H5T_IEEE_F64LE, filespace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    
    	//write data to dataset
    	H5Dwrite (dataset, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT, c->particles);

    	//close dataset, so the next iteration can create it anew with another name
    	H5Dclose (dataset);
    }    
    //close dataspaces   
    H5Sclose (memspace);	//TODO: maybe leave them open and pass them around ?
    H5Sclose (filespace); 
}



//----------------------------------------------------------------------
//close files
static void close_xdmf(struct psc_output_particles_xdmf *out_c) {
	fprintf(out_c->xdmf,"</Grid> \n"); //close temporal grid
	fprintf(out_c->xdmf,"</Domain> \n");
	fprintf(out_c->xdmf,"</Xdmf> \n");
	fclose(out_c->xdmf); 
}
static void close_hdf5(struct psc_output_particles_xdmf *out_c) {
    H5Fclose (out_c->hdf5);
}

// ---------------------------------------------------------------------
// this is what the psc actually calls each step
static void
psc_output_particles_xdmf_run(struct psc_output_particles *out, mparticles_base_t *particles) {
	static int pr;
  	if (!pr) {
  		pr = prof_register("c_out_part", 1., 0, 0);
  	}
  	prof_start(pr);  
  	struct psc *psc = ppsc;
  	struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);

  	if (psc->timestep == 0) {
  		out_c->times_written=0;
  		//ensure expected functionality:
        for(int time=out_c->first; time<ppsc->prm.nmax; time+=out_c->step) {
    		out_c->write_particles[time]=true;
    	}
    	//some counting:
    	for(int time=0; time<ppsc->prm.nmax; time++) {
    		if(out_c->write_particles[time]) {
    			out_c->times_written++;
    		}
    	}
    	create_hdf5(out_c);
    	create_xdmf(out_c);
    	int rank=42;
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	if(rank==0){create_master_xdmf(out_c);}	//the compact master contains no temporal information, so it might as well be written at the start
    	
  	}
  	if (out_c->write_particles[ppsc->timestep]) {
    	write_node_hdf5(out_c,particles);	 
    	write_node_xdmf(out_c,particles);
  	}
  	if (psc->timestep == psc->prm.nmax-1) {
   		close_hdf5(out_c);
    	close_xdmf(out_c);
  	}  
	prof_stop(pr); 
}


//----------------------------------------------------------------------
// every function below this point is legacy stuff from nils' original C output.
// I see no reason to touch any of this, it looks dangerous


static void
psc_output_particles_xdmf_create(struct psc_output_particles *out)
{
  struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);
  out_c->io = mrc_io_create(psc_output_particles_comm(out));
  //I'd like to have the outdir from mrc_io as the default value of data_dir
  mrc_io_set_from_options(out_c->io);
  const char *outdir=NULL;
  mrc_obj_get_param_string((struct mrc_obj*)out_c->io, "outdir" ,&outdir);
  psc_output_particles_set_param_string(out, "data_dir" ,outdir);
  initialize_output(out_c);
}

static void
psc_output_particles_xdmf_setup(struct psc_output_particles *out)
{
  struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);

  const char *outdir=NULL;
  mrc_obj_get_param_string((struct mrc_obj*)out, "data_dir" ,&outdir);

  mrc_io_set_name(out_c->io, "particles_out");
  mrc_io_set_param_string(out_c->io, "basename", "particles");
  mrc_io_set_param_string(out_c->io, "outdir",outdir);
  mrc_io_setup(out_c->io);
  mrc_io_view(out_c->io);
}

static void
psc_output_particles_xdmf_destroy(struct psc_output_particles *out)
{
  struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);
  mrc_io_destroy(out_c->io);
}

static void
psc_output_particles_xdmf_set_from_options(struct psc_output_particles *out)
{
  struct psc_output_particles_xdmf *out_c = to_psc_output_particles_xdmf(out);
  mrc_io_set_from_options(out_c->io);
}

#define VAR(x) (void *)offsetof(struct psc_output_particles_xdmf, x)

static struct param psc_output_particles_xdmf_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },//This default value will NOT be used. The defaul value of the parameter outdir of mrc_io will be used instead.
  { "filter_function"   , VAR(filter_func) , PARAM_PTR(NULL) },
  { "particle_first"       , VAR(first)         , PARAM_INT(0)            },
  { "particle_step"        , VAR(step)          , PARAM_INT(10)           },
  { "write_x"        , VAR(write_variable[0])          , PARAM_BOOL(true)           },
  { "write_y"        , VAR(write_variable[1])          , PARAM_BOOL(true)           },
  { "write_z"        , VAR(write_variable[2])          , PARAM_BOOL(true)           },
  { "write_px"       , VAR(write_variable[3])          , PARAM_BOOL(true)           },
  { "write_py"       , VAR(write_variable[4])          , PARAM_BOOL(true)           },
  { "write_pz"       , VAR(write_variable[5])          , PARAM_BOOL(true)           },
  { "write_charge"   , VAR(write_variable[6])          , PARAM_BOOL(true)           },
  { "write_mass"     , VAR(write_variable[7])          , PARAM_BOOL(true)           },
  { "write_weight"   , VAR(write_variable[8])          , PARAM_BOOL(true)           },
  
  {},
};
#undef VAR


// ======================================================================
// psc_output_particles: subclass "c"

struct psc_output_particles_ops psc_output_particles_xdmf_compact_ops = {
  .name                  = "xdmf_compact",
  .size                  = sizeof(struct psc_output_particles_xdmf),
  .param_descr           = psc_output_particles_xdmf_descr,
  .create                = psc_output_particles_xdmf_create,
  .setup                 = psc_output_particles_xdmf_setup,
  .set_from_options      = psc_output_particles_xdmf_set_from_options,
  .destroy               = psc_output_particles_xdmf_destroy,
//  .write                 = psc_output_particles_c_write,
//  .read                  = psc_output_particles_c_read,
//  .view                  = psc_output_particles_c_view,
  .run                   = psc_output_particles_xdmf_run,
};
