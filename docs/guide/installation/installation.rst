
Installation
************

Getting the code from github
============================

PSC is available from github at
https:://github.com/psc-code/psc.git. To use it, make a local clone
using ``git``:

.. code-block:: sh

   [kai@mbpro ~]$ git clone https://github.com/psc-code/psc.git

Dependencies
============

PSC has the following dependencies:

- (required) `cmake <https://www.cmake.org>`_ is used as the build system

- (required) `MPI
  <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_, the
  Message Passing Interface library used for running the code in parallel

- (included) `googletest <https://github.com/google/googletest>`_ for testing

- ``HDF5``, one option for output

- ``ADIOS2``, another option for ouput and checkpointing

- ``viscid`` is useful for analyzing / visualizing PSC data

.. todo::
   add urls
   
Building
========

PSC uses ``cmake`` as a build system. The build happens in a separate
build directory, e.g. called ``build``:

.. code-block:: sh

   [kai@macbook ~]$ cd psc
   [kai@macbook ~]$ mkdir build
   [kai@macbook ~]$ cd build
   [kai@macbook ~]$ cmake -DCMAKE_BUILD_TYPE=Release ..
   [...]

.. todo::
   add description of cmake options

Hopefully, the cmake step will succeed without error, at which point
we're ready to actually compile the code.

From now on, one should only ever need to build the code using

.. code-block:: sh

   [kai@macbook ~]$ make
   [...]

.. todo::
   psc_harris_xz doesn't compile / work w/o vpic
   
.. todo::
   psc_flatfoil_yz change to small 2d run?

.. todo::
   viscid for viz (incl. mpl3 is_string_view issue)

.. todo::
   add flatfoil_yz notebook etc
   
   
.. todo::
   describe running the tests
