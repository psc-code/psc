
Installation
************

Introduction
============

PSC is available from github at https:://github.com/psc-code/psc.git
. It uses a cmake-based build, so building the code should be
relatively straightforward.

Dependencies
============

PSC has the following dependencies:

- (required) `cmake <https://www.cmake.org>`_ is used as the build system

- (required) `MPI
  <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_, the
  Message Passing Interface library used for running the code in parallel

- (included) `googletest <https://github.com/google/googletest>`_ for testing

- (equired) `HDF5 <https://www.hdfgroup.org/>`_, currently still the only real option
  for regular field / particle output

- (optional) `ADIOS2 <https://github.com/ornladios/ADIOS2>`_, currently the only
  option for checkpoint I/O

- (optional) `viscid <https://viscid-hub.github.io/Viscid-docs/docs/dev/>`_ is
  useful for analyzing / visualizing PSC data


Build Instructions
==================
  
.. toctree::
   :maxdepth: 2

   generic
   summit
   viscid

   

  
