
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda.cxx \
	psc_particles_cuda.cxx \
	psc_push_particles_cuda.cxx \
	psc_push_fields_cuda.cxx \
	psc_marder_cuda.cxx \
	psc_bnd_cuda.cxx \
	psc_bnd_particles_cuda.cxx \
	psc_bnd_fields_cuda.cxx \
	psc_output_fields_item_flds_cuda.cxx \
	psc_output_fields_item_moments_1st_nc_cuda.cxx \
	psc_heating_cuda.cxx \
	psc_inject_cuda.cxx \
	b40c/kernel_utils.h \
	b40c/radixsort_kernel_common.h \
	b40c/radixsort_key_conversion.h \
	b40c/radixsort_reduction_kernel.h \
	b40c/radixsort_scanscatter_kernel.h \
	b40c/radixsort_scanscatter_kernel3.h \
	b40c/radixsort_scanscatter_kernel4.h \
	b40c/radixsort_spine_kernel.h \
	b40c/vector_types.h

libsubdir_la_LIBADD = \
	cuda_base.lo \
	cuda_mparticles_iface.lo \
	cuda_mparticles.lo \
	cuda_mparticles_base.lo \
	cuda_mparticles_bnd.lo \
	cuda_mparticles_bnd_scan.lo \
	cuda_mparticles_bnd_spine.lo \
	cuda_mfields_iface.lo \
	cuda_mfields.lo \
	cuda_mfields_bnd.lo \
	cuda_push_mprts_yz.lo \
	cuda_push_fields_yz.lo \
	cuda_bnd.lo \
	cuda_heating.lo \
	cuda_moments_1st_nc.lo

check_PROGRAMS = \
	tests/test_cuda_accel \
	tests/test_cuda_mfields_1 \
	tests/test_cuda_mparticles_1

#	tests/test_cuda_mparticles
#	tests/test_cuda_heating

#tests_test_cuda_mparticles_SOURCES =
#tests_test_cuda_mparticles_LDADD = tests/test_cuda_mparticles.o cuda_mparticles.lo \
#  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo
#nodist_EXTRA_tests_test_cuda_mparticles_SOURCES = dummy.cxx

tests_test_cuda_accel_SOURCES = 
tests_test_cuda_accel_LDADD = tests/test_cuda_accel.lo \
  cuda_mfields.lo \
  ../../libmrc/src/mrc_json.lo ../../libmrc/src/mrc_json_parser.lo \
  ../../libmrc/src/json-builder.lo ../../libmrc/src/json.lo
nodist_EXTRA_tests_test_cuda_accel_SOURCES = dummy.cxx

tests_test_cuda_mparticles_1_SOURCES = 
tests_test_cuda_mparticles_1_LDADD = tests/test_cuda_mparticles_1.lo \
  cuda_mparticles.lo cuda_mparticles_base.lo \
  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo \
  ../../libmrc/src/mrc_json.lo ../../libmrc/src/mrc_json_parser.lo \
  ../../libmrc/src/json-builder.lo ../../libmrc/src/json.lo
nodist_EXTRA_tests_test_cuda_mparticles_1_SOURCES = dummy.cxx

tests_test_cuda_mfields_1_SOURCES = 
tests_test_cuda_mfields_1_LDADD = tests/test_cuda_mfields_1.lo \
  cuda_mfields.lo \
  cuda_push_fields_yz.lo \
  ../../libmrc/src/mrc_json.lo ../../libmrc/src/mrc_json_parser.lo \
  ../../libmrc/src/json-builder.lo ../../libmrc/src/json.lo
nodist_EXTRA_tests_test_cuda_mfields_1_SOURCES = dummy.cxx

#tests_test_cuda_heating_SOURCES = 
#tests_test_cuda_heating_LDADD = tests/test_cuda_heating.o cuda_mparticles.lo cuda_heating.lo \
#  cuda_mparticles_bnd.lo cuda_mparticles_bnd_scan.lo cuda_mparticles_bnd_spine.lo
#nodist_EXTRA_tests_test_cuda_heating_SOURCES = dummy.cxx

CUDACC_FLAGS = -std=c++11 -g --ptxas-options=-v --maxrregcount=32 --use_fast_math $(filter-out -fopenmp -msse2 -msse4.1 -Wall -ffast-math,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)) $(CUDAMPIFLAGS) $(CUDACFLAGS)

CUDA_DEPS = \
        cuda_mparticles.h \
        cuda_mfields.h \
	cuda_bits.h \
	cuda_iface.h  \
	psc_particles_cuda.h 

cuda_moments_1st_nc.lo: cuda_mparticles_const.h
cuda_push_mprts_yz.lo: cuda_mparticles_const.h
cuda_heating.lo: cuda_mparticles_const.h
cuda_transfer.lo: cuda_iface_bnd.h
cuda_bnd.lo: cuda_iface_bnd.h

#FIXME, this is just really too hacky...

%.lo: %.cu $(CUDA_DEPS)
	@echo "Tricking libtool"
	TMPFILE=${<:.cu=.cu.c}
	echo "void lotrickdummy() {}" > ${<:.cu=.cu.c}
	depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.o$$||'`; \
	$(LTCOMPILE) -MT $@ -MD -MP -MF $$depbase.Tpo -c -o $@ ${<:.cu=.cu.c}
	rm ${@:.lo=.o} ${<:.cu=.cu.c}
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

%.o: %.cu 
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o $@

CLEANFILES = *.fatbin.c *.cudafe1.c *.cudafe2.c *.stub.c *.cubin *.i *.ii *.ptx *.hash *.fatbin *.cudafe1.cpp *.cudafe1.gpu *.cudafe2.gpu *.cu.cpp

endif USE_CUDA

