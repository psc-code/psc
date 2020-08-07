
Installation
************

Introduction
============

PSC is available from github at https:://github.com/psc-code/psc.git
. The recommended build method uses Spack package manager.  Spack is a package manager that targets HPC use cases.  Installing with spack will automate the installation and linking of dependencies. The Legacy build method uses cmake, but the dependencies must be installed manually.

Dependencies
============

PSC has the following dependencies:

- (required) `cmake <https://www.cmake.org>`_ is used as the build system

- (required) `MPI
  <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_, the
  Message Passing Interface library used for running the code in parallel

- (included) `googletest <https://github.com/google/googletest>`_ for testing

- (required) `HDF5 <https://www.hdfgroup.org/>`_, currently still the only real option
  for regular field / particle output

- (optional) `ADIOS2 <https://github.com/ornladios/ADIOS2>`_, currently the only
  option for checkpoint I/O

- (optional) `viscid <https://viscid-hub.github.io/Viscid-docs/docs/dev/>`_ is
  useful for analyzing / visualizing PSC data

- (optional) `rmm <https://github.com/rapidsai/rmm>`_, Performance optimization, provides
  custom cuda allocator to reduce cudaMalloc/cudaFree overhead.

- (optional) `nvtx_pmpi <https://github.com/NVIDIA/cuda-profiler/tree/master/nvtx_pmpi_wrappers>`_ MPI
  hooks so that all MPI calls appear on nsight systems profiler


Spack Build Instructions
========================
- Clone Spack
.. code-block:: sh

   $ git clone -b releases/v0.14 https://github.com/spack/spack.git
- Enable shell support for Spack.

.. code-block:: sh

  # For bash/zsh users
  $ export SPACK_ROOT=/path/to/spack
  $ . $SPACK_ROOT/share/spack/setup-env.sh
- Add the PSC package repo to Spack
.. code-block:: sh

  $ spack repo add path/to/psc/spack/psc
  ==> Added repo with namespace 'psc'.
- And finally,
.. code-block:: sh

   $ spack install psc


Legacy Build Instructions
=========================
                                           
.. toctree::
   :maxdepth: 2

   generic
   summit
   trillian
   viscid


 

  
