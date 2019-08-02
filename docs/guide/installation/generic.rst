
Generic Build Instructions
**************************

Getting the code from github
============================

PSC is available from github at
https:://github.com/psc-code/psc.git. To use it, make a local clone
using ``git``:

.. code-block:: sh

   [kai@mbpro ~]$ git clone https://github.com/psc-code/psc.git

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
we're ready to actually compile the code:

.. code-block:: sh

   [kai@macbook ~]$ make
   [...]

From now on, after making changes, one should only ever need to
rebuild the code using the above ``make`` command.

Running the tests
=================

The PSC code base includes a bunch of unit tests, though coverage is still far from complete. To run these tests, use ``ctest``:

.. code-block:: sh

   [kaig1@login3 build-summit]$ ctest .
   [...]
   
   100% tests passed, 0 tests failed out of 233

   Total Test time (real) = 120.68 sec

