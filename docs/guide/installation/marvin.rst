
Installing on Marvin
*****************************************

This is a temporary solution to installing PSC on Marvin. A better solution would use a spack environment.

Starting from scratch
=====================

First we need git and spack:

.. code-block:: sh

   $ module load git-gcc8/2.24.0
   $ git clone https://github.com/spack/spack.git

To enable spack commands, run

.. code-block:: sh

   $ . ${HOME}/spack/share/spack/setup-env.sh

or just add the line to ``.bash_profile``.

To clone PSC, one may need to set up an RSA key (via ``ssh-keygen``; hit "enter" on all prompts). Copy the full contents of ``id_rsa.pub`` (ie, the output of ``cat`` below) to a new RSA key on GitHub.

.. code-block:: sh

   $ ssh-keygen
   $ cat .ssh/id_rsa.pub

It should now be possible to do 

.. code-block:: sh

   $ git clone git@github.com:psc-code/psc.git
   $ git clone git@github.com:psc-code/psc-infrastructure.git

Continuing with spack...

.. code-block:: sh

   $ spack repo add ${HOME}/psc-infrastructure/spack/psc

At this point, copy the files in ``psc-infrastructure/TODO`` to spack (TODO: currently, one needs to copy the contents of .spack from another user):

.. code-block:: sh

   $ cp -a psc-infrastructure/TODO/* .spack/

Finally, to complete installation of PSC, run

.. code-block:: sh

   $ spack install psc

(this will take a while).

Diagnostics
===========

See the `spack docs <https://spack.readthedocs.io/en/latest/>`_ for how to use spack. Useful commands include ``spack spec``, ``spack compiler list``, and ``spack find``. The ``spack --version`` used here was 0.16.1.


Building
========

In ``psc``, run 

.. code-block:: sh

   $ cmake -B build

Make sure it builds with gcc 10.2. If it fails, most or all of the following may need to be run (or installed, e.g. ``spack install git``). This may not be an exhaustive list. If CMake fails to find anything, identify the relevant package with ``spack find`` and load them as below.

.. code-block:: sh

   $ spack load cmake
   $ spack load git
   $ spack load openmpi
   $ spack load adios2
   $ spack load hdf5
   $ spack load gcc