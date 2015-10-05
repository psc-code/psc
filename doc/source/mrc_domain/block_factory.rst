==============
Block Factory
==============

Class Overview
===============

Block factories are helper classes attached to :c:type:`mrc_domain`
objects of type :c:type:`mb <mrc_domain_mb>`. They are responsible for
carving the global domain into separate blocks and creating the
mappings between those blocks which allow boundary exchanges.

While each block has only three coordinate axes the total number of
computational coordinates in a multi-block domain can be arbitrarily
large. To account for this, block factories are also responsible for
uniquely naming the :c:type:`generators <mrc_crds_gen>` for each
computational coordinate and assigning them to the appropriate axes in
each block.

Because multi-block domains are generally non-Cartesian the idea of a
global number of grid points in each of three dimensions doesn't
apply. As such, each block factory has options to set the number of
grid points along each computational coordinate.

.. warning:: Limitations of the :c:type:`field <mrc_fld>` storage and
	     the domain decomposition require that every block have
	     the same number of grid points in each local coordinate
	     axis, i.e. the logical dimensions of every block must be
	     identical! Block factories should probably check that
	     this is the case...

.. _block_factory_types:

Types
-----

Block diagrams indicate how the domain is split into blocks. The
mapping table shows which faces are connected for each block.

.. c:type:: mrc_block_factory_simple2d

   *Name*: "simple2d"

   *Block Diagram*:
     .. image:: images/simple2d_block_map.png

   *Mapping Table*:
     +-----+------+-----+
     |Block|Face  |Maps |
     |     |      |to   |
     +-----+------+-----+
     | 0   |Left  |N/A  |
     +-----+------+-----+
     |     |Right |N/A  |
     +-----+------+-----+
     |     |Bottom|N/A  |
     +-----+------+-----+
     |     |Top   |N/A  |
     +-----+------+-----+

     .. note:: If a dimension is specified as periodic the appropriate
	       mappings between faces on block 0 will be applied.
   
   *Coord Gen Names*:

     * "coord_gen_x"

     * "coord_gen_y"


   :param param_float [x,y,z]b: Beginning of the domain in [x,y,z]

   :param param_float [x,y,z]e: End of the domain in [x,y,z]

   :param param_int3 m: Number of point in [x,y,z]. :math:`z` direction is ignored and set to 1.

.. c:type:: mrc_block_factory_simple3d

   *Name*: "simple3d"

   *Block Diagram*:
     .. image:: images/simple3d_block_map.png

   *Mapping Table*:
     +-----+------+-----+
     |Block|Face  |Maps |
     |     |      |to   |
     +-----+------+-----+
     | 0   |Left  |N/A  |
     +-----+------+-----+
     |     |Right |N/A  |
     +-----+------+-----+
     |     |Bottom|N/A  |
     +-----+------+-----+
     |     |Top   |N/A  |
     +-----+------+-----+
     |     |Front |N/A  |
     +-----+------+-----+
     |     |Back  |N/A  |
     +-----+------+-----+



     .. note:: If a dimension is specified as periodic the appropriate
	       mappings between faces on block 0 will be applied.


   *Coord Gen Names*:

     * "coord_gen_x"

     * "coord_gen_y"

     * "coord_gen_z"


   :param param_float [x,y]b: Beginning of the domain in [x,y,z]

   :param param_float [x,y]e: End of the domain in [x,y,z]

   :param param_int3 m: Number of point in [x,y,z].


.. c:type:: mrc_block_factory_cylindrical

   *Name*: "cylindrical"

   *Block Diagram*:
     .. image:: images/cylindrical_block_map.png


   *Mapping Table*:
     +----------+----------+-----------+
     |Block     |Face      |Maps To    |
     +----------+----------+-----------+
     |0         |Left      |Block 1    |
     |          |          |Left       |
     |          |          |           |
     +----------+----------+-----------+
     |          |Right     |N/A        |
     +----------+----------+-----------+
     |          |Bottom    |Block 1    |
     |          |          |Top        |
     +----------+----------+-----------+
     |          |Top       |Block 1    |
     |          |          |Bottom     |
     |          |          |           |
     +----------+----------+-----------+
     |1         |Left      |Block 0    |
     |          |          |Left       |
     +----------+----------+-----------+
     |          |Right     |N/A        |
     +----------+----------+-----------+
     |          |Bottom    |Block 0    |
     |          |          |Top        |
     +----------+----------+-----------+
     |          |Top       |Block 0    |
     |          |          |Bottom     |
     +----------+----------+-----------+
     
     .. note:: The left faces for both blocks have a special boundary
	       type `BTYPE_SP` which indicates they lie near a
	       coordinate singularity and must be handled with special care.

   *Coord Gen Names*:

     * "coord_gen_r"

     * "coord_gen_th"


   :param param_float rb: Beginning of the domain in :math:`r`
      
   :param param_float re: End of the domain in :math:`r` (radius of the cylinder)
			  
   :param param_float z[b,e]: Beginning and end of the domain in :math:`z`.

   :param param_int m_r: Number of points in :math:`r`
    
   :param param_int m_th: Number of points in :math:`\theta`

   .. warning:: Setting `rb` to zero could result in NaNs. Some small but finite value is recommended.




Relationships
-------------
* :c:type:`mrc_domain`

  Block Factories exist only to create the block layout for
  :c:type:`"Multi-Block" <mrc_domain_mb>` domains. "mb" domains create,
  call, and destroy block factories. The factory feeds back in a list of
  blocks which the domain uses. The block factory object can be recovered
  from the domain via the :c:func:`get_trafo` method.

* :c:type:`mrc_crds`

  Multi-block domains require coordinates to be generated in a
  consistent manner across all pieces of the domain. Thus block
  factories must store enough :c:type:`mrc_crds_gen` objects to cover
  the domain's different coordinate axes, as well as assign them to
  the appropriate axes in each block. These will be processed by the
  :c:type:`mrc_crds_mb` sub-type to form the domain.

  .. note:: To make it explicit that we're dealing with multi-block
	    domain coordinates (of which there can be more than three)
	    I've adopted the convention of using `coord` instead of
	    `crds`. We'll see if that sticks...

* :c:type:`mrc_trafo`

  Informal connection. Trafos may suggest a block factory to the
  domain. 

User Interface
================

.. c:type:: struct mrc_block_factory

   Block factories are responsible for constructing the domain. The
   concepts of global dimensions and length don't mean much in
   multi-block arrangements. Each block factory will provide some
   public parameter interface for setting the number of grid points
   and coordinate boundaries. Additionally, each block factory will
   define a set of `coord_gen` names which can
   be used to apply nonuniform grids. These are :c:type:`mrc_crds_gen`
   objects except that the names have been changed to indicate that
   they are block coordinates. See the :ref:`types list <block_factory_types>` 
   for what options each type can support.



Writing A Subclass
===================

Required Elements
------------------

.. c:function:: void run(struct mrc_block_factory *fac, struct mrc_domain *domain)

   Assemble the list of blocks for the domain.

   :param struct mrc_block_factory* fac: This block factory
     
   :param struct mrc_domain* domain: The domain where blocks should be
				     assembled. The domain is assured
				     to be multi-block (an assert in
				     the super-class will trip if it's
				     not).

The run function is responsible for allocating space in the `mb` domain
sub-context at the member `mb->mb_blocks`, setting the number of blocks
at `mb->nr_blocks`, then constructing the mapping for each block.

.. warning:: All blocks must have the same number of grid points in
	     the logical dimensions (`0,1,2`), due to the
	     limitations of :c:type:`mrc_fld` storage and domain decomposition.
	     This *does not* mean that mappings between the logical
	     axis and computational coordinate must be the same on
	     each block.


* **Coordinate Generator Objects**: Each computational coordinate in the
  domain must have a :c:type:`mrc_crds_gen` object which is assigned to
  the appropriate block local axis during the run function. It is
  recommended, but not required, that these objects be registered as
  member objects and have sensible public names so they can be
  manipulated from the command line.

Optional Elements
------------------
Any input parameters which are necessary to initialize the block
mappings and domain, such as numbers of grid points, 


Tutorial
--------------

Constructing block factory sub-types can be rather complicated. This
tutorial will examine, in depth, the construction of the 
:c:type:`cylindrical <mrc_block_factory_cylindrical>` block
factory, which is a relatively simple two block factory with three
global computational coordinates. 


Input Parameter Definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We need to define a subclass context to contain the public
parameters. This should have elements for, at the least, the
coordinate dimensions, bounds, and generation objects. For the
cylindrical domain the bounds on the :math:`\theta` coordinate are
fixed at :math:`[0, 2\pi]`, so we only need bounds in :math:`r` and
:math:`z`. Currently this factory only creates 2D domains and
the number of grid points in :math:`\theta` is split evenly between
the two blocks, so we only need two dimensions variables and will
hard-code the number of :math:`z` grid points to 1. Despite this, we
still need three :c:type:`coord_gen <mrc_crds_gen>` objects, one for
each of the three coordinates. The :c:type:`coordinate object <mrc_crds>` 
will still expect a `coord_gen` for each block local direction, even
if that direction only has one point. These generator objects should
be registered as member objects (`MRC_VAR_OBJ`) in the parameter
definition so they will be created, destroyed, and setup along with the subclass.

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 35-58

Setting the `coord_gen` Names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subclass create function is required to set the coordinate generator
names. As names will be used on the command line to modify grid point
distributions the names should be unique. The following convention is
suggested: `coord_gen_$` where `$` is the name of the computational
coordinate.

.. warning:: A read function which reads member objects *must* exist if the subclass create is
	     used to set the generator names, lest the subclass create
	     try to set the names during read in when the
	     member objects haven't been created. Yes, this is
	     hacky. No, I don't care enough to fix it.

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 60-75


The `run` Function
^^^^^^^^^^^^^^^^^^^^

Housekeeping
"""""""""""""

The run function needs to assign the number of blocks in the
:c:type:`mrc_domain_mb` sub-context and allocate the space to hold
them.

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 82-90
	   
This function should also initialize the block mappings into the
allocated space.

Accessing Domain Parameters
""""""""""""""""""""""""""""

Block factories may assume that all the publicly tunable parameters of
the domain have been set, such as the `mb->bc` boundary condition
parameters. These can be used to set the block mappings. As block
factories are run at the beginning of the `mrc_domain` setup process
one *cannot* assume private members of the domain have been
initialized.



Block Mappings
^^^^^^^^^^^^^^^

Each block contains a lot of information as to where it fits in the
domain.

.. c:type:: struct MB_block

   .. c:member:: int nr_block 
      
      The number of this block (seems redundant but *must* be set)
  
   .. c:member:: int mx[3] 

      Number of points in each block local logical dimension

   .. c:member:: struct MB_face faces[NR_FACES]
    
      Face mappings (6 total)
   
   .. c:member:: struct mrc_crds_gen* coord_gen[3]
   
      Coordinate generators for the block local coordinates.

   .. c:member:: float xl[3]
   
      Lower bounds of this block in coordinate space (passed to the
      `coord_gen` objects)

   .. c:member:: float xh[3]
     
      Upper bounds of this block in coordinate space. (passed to the
      `coord_gen` objects)

For the first block of our cylindrical example, the initialization looks
like this:

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 119-139

And the second block:

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 141-157


Let's go over each of these elements in more depth.

Block Number
"""""""""""""

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 120

Assigning the number of this block is sort of redundant, but it's need
for backwards compatibility with some of the domain systems.

Grid Dimensions
""""""""""""""""

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 121

Each factory should have parameters allowing the grid size in each
computational coordinate to be specified. These are then used to set
the block dimensions in the `.mx` attribute. Note that all blocks must
be the same size even if they don't share computational
coordinates. This restriction comes from the storage of
:c:type:`mrc_fld` data and the :c:type:`mrc_domain_mb` decomposition
methods. It would probably be a good idea to add some checks to each
factory to ensure this is the case.

Face Mappings
""""""""""""""

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 122-133

Each block is a standard 3D cube with faces:

  * FACE_LEFT (:math:`-X`)

  * FACE_RIGHT (:math:`+X`)

  * FACE_BOTTOM (:math:`-Y`)

  * FACE_TOP (:math:`+Y`)

  * FACE_BACK (:math:`-Z`)

  * FACE_FRONT (:math:`+Z`)

Any face not mapped is assumed to be a symmetry direction. In the
example above the `Z` axial dimension of the cylinder is handled via a
helical symmetry elsewhere in the code, so the faces aren't mapped.

Each face requires a `.map` element which aligns the `X`, `Y`, and `Z`
axes on that face with the axes on the target face. For example, the
map of the block 0 bottom face is `{ MB_X, 0, MB_Z }`, which
indicates that the `X` and `Z` axes of the block 0 bottom face align
with the `X` and `Z` axes of the target (block 1 top) face, and the
`Y` direction is normal to the face (ie, the direction in which ghost
points should be exchanged). Map directions may be any of the
following:

  * MB_NONE = 0: Face normal
  
  * MB_X

  * MB_Y

  * MB_Z 

  * MB_XR : X direction but reversed

  * MB_YR : Y direction but reversed

  * MB_ZR : Z direction but reversed

Additionally each face should specify a boundary type
`.btype`. Options are:

  * BTYPE_NONE : The default simple ghost exchange

  * BTYPE_OUTER : Outer boundary, no boundary exchange occurs

  * BTYPE_SP : The boundary is a coordinate singularity. Ghost point
    exchange will occur but additional processing must be performed.

  * BTYPE_SPC : The boundary is a transition between two coordinate
    systems. Ghost point exchange will occur but additional processing
    must be performed.

  * BYTPE_USER : The boundary is unique. Ghost point exchange will
    occur but user generated code will be needed to finalize the boundary.


Coordinate Generation
""""""""""""""""""""""

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 134-138

Each factory sub-type should have, as member objects, enough
:c:type:`coordinate generators <mrc_crds_gen>` to handle its
computational coordinates. The generators should have unique names set
during sub-type creation that would allow command line parsing of
coordinate generation parameters, otherwise there will be issues with
generating non-uniform coordinates. I've adopted the convention of
using the `coord` prefix, instead of the standard `crds`, to indicate
that these code generators are associated with computational
coordinates in a multi-block domain, but that may change. 

These generators must be assigned to the appropriate local axes of each
block in the `.coord_gen` attribute. Additionally, the lower and upper
coordinate bounds of each block are set in the `.xl` and `.xh`
attributes, respectively.

Wrapping it up
^^^^^^^^^^^^^^^^

All that remains after the run function is a standard
:c:type:`mrc_obj_ops` declaration.

.. literalinclude:: ../../../src/mrc_block_factory_cylindrical.c
   :lines: 161-171


Keep in mind that the sub-type must be registered before it can be
selected:

.. literalinclude:: ../../../src/mrc_block_factory.c
   :lines: 89

