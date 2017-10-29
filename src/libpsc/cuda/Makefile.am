
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda.c \
	psc_particles_cuda.c \
	psc_push_particles_cuda.c \
	psc_push_fields_cuda.c \
	psc_marder_cuda.c \
	psc_bnd_cuda.c \
	psc_bnd_particles_cuda.c \
	psc_bnd_fields_cuda.c \
	psc_output_fields_item_flds_cuda.c \
	psc_output_fields_item_moments_1st_nc_cuda.c \
	psc_heating_cuda.c \
	psc_inject_cuda.c \
	b40c/kernel_utils.h \
	b40c/radixsort_kernel_common.h \
	b40c/radixsort_key_conversion.h \
	b40c/radixsort_reduction_kernel.h \
	b40c/radixsort_scanscatter_kernel.h \
	b40c/radixsort_scanscatter_kernel3.h \
	b40c/radixsort_scanscatter_kernel4.h \
	b40c/radixsort_spine_kernel.h \
	b40c/vector_types.h \
	cuda_common.h \
	cuda_sort2.h \
	cuda_sort2_spine.h \
	particles_cuda.h \
	psc_bnd_cuda.h \
	psc_bnd_cuda_fields.h \
	psc_cuda.h

libsubdir_la_LIBADD = \
	cuda_base.lo \
	cuda_mparticles.lo \
	cuda_mparticles_bnd.lo \
	cuda_mparticles_bnd_scan.lo \
	cuda_mparticles_bnd_spine.lo \
	cuda_mfields.lo \
	cuda_transfer.lo \
	cuda_push_mprts_yz.lo \
	cuda_push_fields_yz.lo \
	cuda_bnd.lo \
	cuda_spine.lo \
	cuda_exclusive_scan_2.lo \
	cuda_exchange_particles.lo \
	cuda_heating.lo \
	cuda_moments_1st_nc.lo

#	cuda_find_indices.lo \
	cuda_sort.lo \
	cuda_sort2.lo

check_PROGRAMS = \
	tests/test_cuda_mfields_1 \
	tests/test_cuda_mparticles_1

#	tests/test_cuda_mparticles
#	tests/test_cuda_heating

tests_test_cuda_mparticles_SOURCES =
tests_test_cuda_mparticles_LDADD = tests/test_cuda_mparticles.o cuda_mparticles.lo \
  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo
nodist_EXTRA_tests_test_cuda_mparticles_SOURCES = dummy.cxx

tests_test_cuda_mparticles_1_LDADD = cuda_mparticles.lo \
  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo \
  ../../libmrc/src/mrc_json.lo ../../libmrc/src/mrc_json_parser.lo \
  ../../libmrc/src/json-builder.lo ../../libmrc/src/json.lo
nodist_EXTRA_tests_test_cuda_mparticles_1_SOURCES = dummy.cxx

tests_test_cuda_mfields_1_LDADD = cuda_mfields.lo \
  cuda_push_fields_yz.lo \
  ../../libmrc/src/mrc_json.lo ../../libmrc/src/mrc_json_parser.lo \
  ../../libmrc/src/json-builder.lo ../../libmrc/src/json.lo
nodist_EXTRA_tests_test_cuda_mfields_1_SOURCES = dummy.cxx

tests_test_cuda_heating_SOURCES = 
tests_test_cuda_heating_LDADD = tests/test_cuda_heating.o cuda_mparticles.lo cuda_heating.lo \
  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo
nodist_EXTRA_tests_test_cuda_heating_SOURCES = dummy.cxx

CUDACC_FLAGS = -g --ptxas-options=-v --maxrregcount=32 --use_fast_math $(filter-out -fopenmp -msse2 -msse4.1 -Wall -ffast-math,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)) $(CUDAMPIFLAGS) $(CUDACFLAGS)

CUDA_DEPS = \
	cuda_iface.h \
        cuda_mparticles.h \
        cuda_mfields.h \
	cuda_bits.h \
        psc_cuda.h \
        constants.c

#        psc_particles_cuda.h FIXME

cuda_spine.lo: b40c/radixsort_reduction_kernel.h b40c/radixsort_scanscatter_kernel.h cuda_sort2_spine.h
cuda_sort2.lo cuda_exchange_particles.lo cuda_sort2_spine.lo: psc_cuda.h cuda_sort2.h
cuda_exchange_particles.lo: psc_bnd_cuda.h ../../include/psc_particles_cuda.h cuda_mparticles_const.h
cuda_exclusive_scan_2.lo: ../../include/psc_particles_cuda.h
cuda_transfer.lo: ../../include/psc_particles_cuda.h
cuda_spine.lo: ../../include/psc_particles_cuda.h b40c/radixsort_reduction_kernel.h b40c/radixsort_scanscatter_kernel3.h
cuda_moments_1st_nc.lo: cuda_mparticles_const.h
cuda_push_mprts_yz.lo: cuda_mparticles_const.h
cuda_heating.lo: cuda_mparticles_const.h

#FIXME, this is just really too hacky...

%.lo: %.cu $(CUDA_DEPS)
	@echo "Tricking libtool"
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm ${@:.lo=.o}
	rm dummy.c
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

%.o: %.cu 
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o $@

CLEANFILES = *.fatbin.c *.cudafe1.c *.cudafe2.c *.stub.c *.cubin *.i *.ii *.ptx *.hash *.fatbin *.cudafe1.cpp *.cudafe1.gpu *.cudafe2.gpu *.cu.cpp

endif USE_CUDA

