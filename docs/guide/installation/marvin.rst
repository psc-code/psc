
Building and Running on Marvin/Plasma (at UNH)
**********************************************

There are two main approaches using Spack. The first option is using Spack to directly install psc (i.e., ``spack install psc``) after a few preparatory steps, and then later manually loading PSC's immediate dependencies to actually compile (e.g., ``spack load cmake``). The second approach is outlined below and results in a Spack environment in which PSC can be compiled.

Loading Git
===========

.. code-block:: bash
   $ module load git-gcc8

Installing Spack
================

In your home directory:

.. code-block:: bash
   $ git clone https://github.com/spack/spack.git

Note: the version of Spack used here was commit 2d77e44f6f57dcd1529c1ce712acfc6a60f47206, "Pcluster local buildcache". Release 1.4.1 should also be fine. You can swich spack versions by checkout out the desired commit in the spack repo.

Enable ``spack`` commands by running:

.. code-block:: bash
   $ . ~/spack/share/spack/setup-env.sh

It is suggested that you add the above line to your ``.bash_profile``.

Fixing Spack
============

Spack fails to do several things out-of-the-box. First, fix ``~/.spack/cray/compilers.yaml`` by specifying the paths for your compiler. The default is gcc8.2.0, as below. Either delete all other compilers in this yaml, specify the appropriate paths for all of the compilers, or specify ``%gcc@8.2.0`` in your Spack specs.

.. code-block:: yaml
   # file: ~/.spack/cray/compilers.yaml
   - compiler:
      spec: gcc@=8.2.0
      paths:
         cc: /cm/shared/apps/gcc8/8.2.0/bin/gcc
         cxx: /cm/shared/apps/gcc8/8.2.0/bin/g++
         f77: /cm/shared/apps/gcc8/8.2.0/bin/gfortran
         fc: /cm/shared/apps/gcc8/8.2.0/bin/gfortran
      flags: {}
      operating_system: rhel7
      target: any
      modules:
         - PrgEnv-gnu
         - gcc/8.2.0
      environment: {}
      extra_rpaths: []

Second, it is important to use Marvin's built-in slurm. Edit ``~/.spack/packages.yaml`` like so:

.. code-block:: yaml
   # file: ~/.spack/packages.yaml
   packages:
     slurm:
       externals:
         - spec: slurm@18.08.9
           prefix: /cm/shared/apps/slurm/18.08.9

Note the prefix (and version) should agree with the results of ``module show slurm``. Specifically, the prefix needs to be a path to a folder containing both ``include/`` and ``lib64/`` (or symlinks thereto) for slurm.

Setting Up Spack Environment
============================

Make a Spack environment named "psc" (or whatever) and edit the config to match the yaml below. As mentioned above, you may want to add ``%gcc@8.2.0`` (or whatever compiler version you're using) to each of these specs.

.. code-block:: bash
   $ spack env create psc
   $ spack -e psc config edit

.. code-block:: yaml
   # file: ~/spack/var/spack/environments/psc/spack.yaml
   spack:
      specs:
         - "cmake@3.17.0:"
         - gtensor
         - hdf5+hl
         - "adios2@2.4.0:2.8.3"
         - "googletest@1.10.0:"
         - "openmpi+pmi schedulers=slurm"
   view: true
   concretizer:
      unify: true

Then concretize and install the packages (this takes a while). There should be a single warning after the install about skipping slurm, an external package.

.. code-block:: bash
   $ spack -e psc concretize
   $ spack -e psc install

Cloning PSC
===========

You may need to set up an RSA key (via ``ssh-keygen``; hit "enter" on all prompts). Copy the full contents of ``id_rsa.pub`` (ie, the output of ``cat`` below) to a new RSA key on GitHub.

.. code-block:: bash
   $ ssh-keygen
   $ cat ~/.ssh/id_rsa.pub

Clone PSC wherever you want it:

.. code-block:: bash
   $ git clone git@github.com:psc-code/psc.git

Compiling PSC
=============

Before compiling the code, you must first activate the spack environment:

.. code-block:: bash
   $ spack env activate psc

This command is slow and may interfere with subsequent slurm commands, so it is *not* recommended that you put this in your bash profile. Just activate it when working on the code.

To actually compile, go into the cloned repository and run

.. code-block:: bash
   $ cmake -B build
   $ cd build
   $ make

Running PSC
===========

NOTE: slurm commands may not work while the spack environment is active. Deactivate it first.

Use slurm to run PSC. For example, run 

.. code-block:: bash
   $ sbatch run.sh

where ``run.sh`` looks like

.. code-block:: bash
   #!/bin/bash
   #SBATCH --nodes=4
   #SBATCH --ntasks-per-node=64
   #SBATCH --cpus-per-task=1
   #SBATCH --time=24:00:00

   cp path/to/psc/bits/adios2cfg.xml adios2cfg.xml
   srun path/to/psc/build/src/executable