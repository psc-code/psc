
Installing on Marvin (at UNH)
*****************************************

Starting from scratch
=====================

First we need git and spack:

.. code-block:: sh

   $ module load git-gcc8/2.24.0
   $ git clone https://github.com/spack/spack.git

To enable spack commands, run

.. code-block:: sh

   $ . ${HOME}/spack/share/spack/setup-env.sh

or just edit ``.bashrc`` as appropriate.

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

At this point, copy the files in ``psc-infrastructure/TODO`` to spack (TODO or symlink them? Is that easier?):

.. code-block:: sh

   $ cp -a psc-infrastructure/TODO/* .spack/

Finally, to complete installation of PSC, run

.. code-block:: sh

   $ spack install psc

(this will take a while).

Diagnostics
===========

To diagnose problems with spack, one may use

.. code-block:: sh

   $ spack spec psc

to see a list of psc's dependencies;

.. code-block:: sh

   $ spack --version

to see the version (mine was 0.16.1);

.. code-block:: sh

   $ spack --help

to see options;

.. code-block:: sh

   $ spack compiler list

to see available compilers.