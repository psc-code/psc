
Building and Running on Trillian (at UNH)
*****************************************

Setting up the environment
==========================

I have a little script that I source after logging in:

.. code-block:: sh

   [kaig1@login3 ~]$ cat setup-psc-env-trillian.sh
   PATH=/home/space/kai/bin:$PATH
   module load gcc/7.3.0

   [kaig1@login3 ~]$ source setup-psc-env-trillian.sh

It might be even better to do these things from your ``.bashrc``, so they'll be available automatically every time you log in.

Getting the code from github
============================

I keep my codes in the ``~/src`` directory, so I'm using this in the following example.

PSC is available from github at
https:://github.com/psc-code/psc.git. To use it, first make a fork on
github into your own account. (Click the ``Fork`` button in the upper
right corner). You then should make a local clone
using ``git``:

.. todo::
   Unify fork / clone instructions.

.. code-block:: sh

   [kai@trillian ~]$ cd src
   [kai@trillian src]$ git clone git@github.com:germasch/psc.git
   Cloning into 'psc'...

Building
========

PSC uses ``cmake`` as a build system. The build happens in a separate
build directory, e.g. called ``build-trillian``:

.. code-block:: sh

   [kaig1@login3 src]$ cd psc
   [kaig1@login3 psc]$ mkdir build-trillian
   [kaig1@login3 psc]$ cd build-trillian

   I create another little script file ``cmake.sh`` in the ``build-trillian/`` directory, so that I know how I invoked ``cmake`` if I need to do it again in the future:

.. code-block:: sh

   #! /bin/bash

   cmake \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-Wno-undefined-var-template" \
   ..

Now, invoke cmake:

.. code-block:: sh

   [kaig1@trillian build-trillian]$ ./cmake.sh
   -- The C compiler identification is GNU 7.3.0
   -- The CXX compiler identification is GNU 7.3.0
   [...]

Hopefully, the cmake step will succeed without error, at which point
we're ready to actually compile the code:

.. code-block:: sh

  [kaig1@trillian build-trillian]$ make		
  [...]

From now on, after making changes, one should only ever need to
rebuild the code using the above ``make`` command.

Running the tests
=================

Running the tests on trillian is kinda non-trivial, since they are
using MPI, so they require to be run with `aprun`, but that only works
if you're inside of a batch job.

Running a job
=============

Here is a job script ``harris.sh`` to run the small sample 2-d flatfoil case on Trillian:

.. code-block:: sh

   #! /bin/bash
   #PBS -l nodes=1:ppn=32
   #PBS -l walltime=00:10:00
    
   DIR=~/scratch/harris/harris_001
   mkdir -p $DIR
   cd $DIR

   cp ~/src/psc/src/psc_harris_xz.cxx .

   aprun -n 4 ~/src/psc/build-trillian/src/psc_harris_xz \
     2>&1 | tee log

Submit as usual:

.. code-block:: sh

   [kaig1@trillian build-trillian]$ qsub harris.sh
   [kai@trillian build-trillian]$ qstat -u kai

   sdb:
                                                                Req'd  Req'd   Elap
   Job ID          Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
   --------------- -------- -------- ---------- ------ --- --- ------ ----- - -----
   57623.sdb       kai      workq    run.sh      22952   1  32    --  00:10 R 00:04
