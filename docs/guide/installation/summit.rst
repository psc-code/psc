
Building and Running on Summit
******************************

Setting up the environment
==========================

I have a little script that I source after logging in:

.. code-block:: sh

   [kaig1@login3 ~]$ cat setup-psc-env-summit.sh
   module load cmake/3.14.2
   module load gcc/7.4.0
   module load hdf5/1.10.3
   module load cuda/10.1.168

   [kaig1@login3 ~]$ source setup-psc-env-summit.sh

It might be even better to load those modules from your ``.bashrc``, so they'll be available automatically every time you log in.

Getting the code from github
============================

I keep my codes in the ``~/src`` directory, so I'm using this in the following example.

PSC is available from github at
https:://github.com/psc-code/psc.git. To use it, make a local clone
using ``git``:

.. code-block:: sh

   [kaig1@login3 ~]$ cd src
   [kaig1@login3 src]$ git clone https://github.com/psc-code/psc.git
   Cloning into 'psc'...
   remote: Enumerating objects: 1245, done.
   remote: Counting objects: 100% (1245/1245), done.
   remote: Compressing objects: 100% (387/387), done.
   remote: Total 102767 (delta 885), reused 1111 (delta 838), pack-reused 101522
   Receiving objects: 100% (102767/102767), 39.34 MiB | 24.27 MiB/s, done.
   Resolving deltas: 100% (84389/84389), done.
   Checking out files: 100% (1467/1467), done.

Building
========

PSC uses ``cmake`` as a build system. The build happens in a separate
build directory, e.g. called ``build-summit``:

.. code-block:: sh

   [kaig1@login3 src]$ cd psc
   [kaig1@login3 psc]$ mkdir build-summit
   [kaig1@login3 psc]$ cd build-summit

I create another little script file ``cmake.sh`` in the ``build-summit/`` directory, so that I know how I invoked ``cmake`` if I need to do it again in the future:

.. code-block:: sh

   #! /bin/bash

   cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx \
   ..

Now, invoke cmake:

.. code-block:: sh

   [kaig1@login3 build-summit]$ ./cmake.sh
   -- The C compiler identification is GNU 7.4.0
   -- The CXX compiler identification is GNU 7.4.0
   [...]

Hopefully, the cmake step will succeed without error, at which point
we're ready to actually compile the code:

.. code-block:: sh

  [kaig1@login3 build-summit]$ make		
  [...]

From now on, after making changes, one should only ever need to
rebuild the code using the above ``make`` command.

.. note::

   On Summit, building the tests creates a bunch of annoying warnings
   since cmake runs those executables to discover the tests, and on Summit,
   running those executables gives warnings. The following helps to quiet those
   down somewhat (but don't use this when actually running the code with mpirun).

   .. code-block:: sh

      [kaig1@login3 build-summit]$ export OMPI_MCA_btl=tcp,self

   

Running the tests
=================

The PSC code base includes a bunch of unit tests, though coverage is still far from complete. To run these tests, use ``ctest``:

.. code-block:: sh

   [kaig1@login3 build-summit]$ ctest .
   [...]
   
   100% tests passed, 0 tests failed out of 233

   Total Test time (real) = 120.68 sec

Running a job
=============

Here is a job script ``flatfoil.sh`` to run the small sample 2-d flatfoil case on Summit:

.. code-block:: sh

   #! /bin/bash
   #BSUB -P AST147
   #BSUB -W 00:10
   #BSUB -nnodes 1
   #BSUB -J flatfoil_summit004

   DIR=$PROJWORK/ast147/kaig1/flatfoil-summit004
   mkdir -p $DIR
   cd $DIR

   jsrun -n 4 -a 1 -c 1 -g 1 ~/src/psc/build-summit-gpu/src/psc_flatfoil_yz

Submit as usual:

.. code-block:: sh

   [kaig1@login3 build-summit-gpu]$ bsub flatfoil.sh
   Job <523811> is submitted to default queue <batch>.
   [kaig1@login3 build-summit-gpu]$ bjobs
   JOBID   USER       STAT   SLOTS    QUEUE       START_TIME    FINISH_TIME   JOB_NAME
   523811  kaig1      PEND      -     batch             -             -       flatfoil_summit004

