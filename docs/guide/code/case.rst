
PSC Cases
*********

A PSC case is a source file that defines a particular simulation run. This corresponds to what otherwise might be refered to as a parameter file or input deck. It is, however, actual C++ code that needs to be compiled when changes are made to it.

The goal is to have PSC cases to be self-contained, ie., everything related to particular run should be in this one file, including plasma parameters, control parameters, boundary conditions, initial conditions, etc. Having cases written in an actual programming language has the advantage of being fully flexible and powerful. There are drawbacks, though:

* In particular, those parameter files are subject to break as PSC evolves. This is currently likely to happen frequently, as major work is still going on with the code.

* One cannot just change a parameter in the case and run it -- the code needs to be recompiled first.

Alternatively, use an InputParams object to load parameters from an input file at runtime. See ``src/include/input_params.hxx``.

Changing / adding a case
========================

The eventual goal is for PSC's API to settle down to a stable state. At that point, the idea is that cases can be maintained externally (e.g., in their own repository), and that an easy way to build those against an installed PSC (library) will be provided. However, as the code is currently far from that goal, one should just adapt cases as part of the main repository. The natural starting point is the existing ``src/psc_flatfoil_yz.cxx`` case. One can either modify that directly, or, if larger changes are to be made copy this file to ``src/psc_my_case.cxx`` and work with that. To add the new ``psc_my_case`` to the build, simply add a line ``add_psc_executable(psc_my_case)`` to ``src/CMakeLists.txt``.

Anatomy of a case
=================

`src/psc_flatfoil_yz.cxx` has been commented to help understand what goes into a case.

.. todo::

   There really should be a description here, too...

   
