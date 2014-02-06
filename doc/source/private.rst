========================
Private(Helper) Classes
========================

Instances of these classes are typically not directly accesible to the
user. They instead exist as private objects associated with one of the
public classes that 'help' the main class by performing some specific
task. In most cases the associated public class will have a parameter
which can be set to change the behavior of the helper class (e.g.
:c:type:`mrc_crds` public objects have options to set which private
:c:type:`mrc_crds_gen` helper objects are used in each dimension.)

For the most part these exist to make the code more module and the
developer's lives a little easier.

.. toctree::
   :maxdepth: 2

   mrc_domain/block_factory
