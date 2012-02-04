PSC Particle Output v0.4
by Simon Jagoda
20.12.2011


Overview:
1. Included output variants
2. Performance comparison
3. General Usage
4. Known issues and stuff left to do



1.	In addition to the C output by Nils, the old/broken fortran output
	and the trusty "none" option, the PSC now offers the following
	choices for particle output:

1.1 xdmf_compact: Each MPI node writes 1 hdf5 container and 1
	corresponding xdmf file, containing the data for ALL output steps.
	These xdmfs are then linked together by file particles_master.xdmf 
	(2*nodes+1 files total). The collection	structure used in the files 
	is spatial-temporal-spatial. Due to a special "feature" in paraview,
	some versions of the software might refuse to open xdmf files with 
	spatial collections of temporal collections (according to google, 
	but works fine with 3.4.0).
	This is intended for people who prefer having few large files over
	thousands of smaller ones.
	
1.2 xdmf_spread: Also using the hdf5/xdmf combination, but with one pair
	of files per output step and node, as well as 1 spatial master file
	per timestep and the temporal master file particles_master.xdmf
	(2*nodes*timesteps+timesteps+1 files total).
	All parts of this output can be used by themselves without the rest
	and the temporal master file is completely rewritten each output step,
	which means you can track progress while the simulation is running.
	Collection structure here is temporal-spatial-spatial, which has no
	known issues. Will usually produce huge amounts of rather small files.
	Suitable for people who like to take a quick look at only a small
	part of their output without having to copy all the data or who
	dislike their filesystem.
	
1.3 custom_binary: Formerly called "c2", dumps binary data using a rather
	unique file structure which will need additional treatment before it
	can be opened by any visualization program known to me. While highly
	impractical to use, it is the most efficient option both in terms of
	output speed and disk usage.
	This one is raw power with no comfort or easy data accessibility
	whatsoever.
	
	 
2. 	Performance 
	Output speed and file sizes have been tested with a simulation on capp2,
	containing 3 kinds of particles each 16000 strong, on 16 patches spread
	over 4 MPI nodes. Speed is averaged over 10 output steps, while size
	is the total amount of space used by the output files at end of simulation.
	
	xdmf_compact:	69.3ms		32.6MB
	xdmf_spread:	83.7ms		36.7MB
	custom_binary:	33.5ms		29.3MB	
	c:				667363.2ms	540.7MB
	
	
3.	General usage

3.1	Choose an output mode: 
	psc_output_particles_set_type(psc->output_particles,"<mode>");	
	in the psc_create part of your case. Working options for <mode>
	are xdmf_compact, xdmf_spread, custom_binary, c and none.
	
3.2 Choose which variables to write: 	
	psc_output_particles_set_param_bool(psc->output_particles, "write_<variable>" ,<true/false>);
	in psc_set_from_options. Available  choices for <variable> are: x, y,
	z, px, py, pz, charge, mass and weight. Default is true for all.
	
3.3 Choose output steps:
	psc_output_particles_set_param_int(psc->output_particles, "particle_step" , <stepsize>);
 	psc_output_particles_set_param_int(psc->output_particles, "particle_first" , <firstoutput>);
 	both in psc_set_from_options. This ensures the expected functionality
 	of output at timesteps with a fixed interval <stepsize>, beginning at
 	<firstoutput>.

3.4 Choose output directory:
	psc_output_particles_set_param_string(psc->output_particles, "data_dir" ,"<path>");
	in psc_set_from_options. Highly recommended for xdmf_spread.
	
	
4. TODO
-	Bug: figure out why paraview writes a "could not read hdf5" error to
	terminal while correctly reading files created by xdmf_compact
-	Bug: get rid of the other paraview warnings as well, if possible
-	Feature: correctly implement output steps on a per timestep level
-	Feature: activate filter functions
-	Misc: get rid of all the compiler warnings	
 	
	
	
	
	
	
