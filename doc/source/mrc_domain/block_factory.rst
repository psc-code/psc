==============
Block Factory
==============

Class Overview
===============

Block factories are helper classes attached to :c:type:`mrc_domain`
objects of type :c:type:`mb <mrc_domain_mb>`. They are responsible for
carving the global domain into seperate blocks and creating the
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


   :param param_float [x,y,z]b: Beginning of the domian in [x,y,z]

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


   :param param_float [x,y]b: Beginning of the domian in [x,y,z]

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
  :c:type:`"Multi-Block" <mrc_domain_mb>` domains. "mb" domains set the type, create,
  and call block factories. The factory feeds back in a list of
  blocks which the domain uses.

* :c:type:`mrc_crds`

  Multi-block domains require coordinates to be generated in a
  consistent manner across all pieces of the domain. Thus block
  factories must store enough :c:type:`mrc_crds_gen` objects to cover
  the domain's differenent coordinate axes, as well as assign them to
  the approriate axes in each block. These will be processed by the
  :c:type:`mrc_crds_mb` subtype to form the domain.

  .. note:: To make it explicit that we're dealing with multi-block
	    domain coordinates (of which there can be more than three)
	    I've adopted the convetion of using `coord` instead of
	    `crds`. We'll see if that sticks...

* :c:type:`mrc_trafo`

  Informal connection. Trafos may suggest a block factory to the
  domain. 

User Interface
================

.. c:type:: struct mrc_block_factory

   Standard block factories are hidden from the user. Type selection
   should be handled purely through the :c:type:`mrc_domain_mb`
   interface. 

   In the future individual subtypes may have command line tunable
   parameters. For this to be supported the domain would need to make
   a `mrc_block_factory_set_from_options` call during its setup
   phase. A minor change.

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
				     the superclass will trip if it's
				     not).

The run function is responsible for allocating space in the `mb` domain
subcontext at the member `mb->mb_blocks`, setting the number of blocks
at `mb->nr_blocks`, then constructing the mapping for each block.

The best way to illustrate the structure of a `run` function is
through this example from the :c:type:`mrc_block_factory_cylindrical`
subtype::

  FIXME: I want to include a source file example here, but I don't
  want to copy paste it. There has to be a better way, I just need to
  learn what it is.



Grid Dimensions
^^^^^^^^^^^^^^^^

Each factory should have parameters allowing the grid size in each
computational coordinate to be specified. These are then used to set
the block dimensions in the `.mx` attribute. Note that all blocks must
be the same size even if they don't share computational
coordinates. This restriction comes from the storage of
:c:type:`mrc_fld` data and the :c:type:`mrc_domain_mb` decomposition
methods. It would probably be a good idea to add some checks to each
factory to ensure this is the case.

Face Mappings
^^^^^^^^^^^^^^

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

  * MB_XR : X direction but reveresed

  * MB_YR : Y direction but reveresed

  * MB_ZR : Z direction but reveresed

Additionally each face should specify a boundary type
`.btype`. Options are:

  * BTYPE_NONE : The default simple ghost exchange

  * BTYPE_OUTER : Outer boundary, no boundary exchange occurs

  * BTYPE_SP : The boundary is a coordinate singularity. Ghost point
    exchange will occur but addtional processing must be performed.

  * BTYPE_SPC : The boundary is a transition between two coordinate
    systems. Ghost point exchange will occur but additional processing
    must be performed.

  * BYTPE_USER : The boundary is unique. Ghost point exchange will
    occur but user generated code will be needed to finalize the boundary.


Coordinate Generation
^^^^^^^^^^^^^^^^^^^^^^

Each factory subtype should have, as member objects, enough
:c:type:`coordinate generators <mrc_crds_gen>` to handle its
computational coordinates. The generators should have unique names set
during subtype creation that would allow command line parsing of
coordinate generation parameters, otherwise there will be issues with
generating non-uniform coordinates. I've adopted the convention of
using the `coord` prefix, instead of the standard `crds`, to indicate
that these code generators are associated with computational
coordinates in a multi-block domain, but that may change. 

These generators must be assigned to the appropriate local axes of each
block in the `.coord_gen` atrribute. Additionally, the lower and upper
coordinate bounds of each block are set in the `.xl` and `.xh`
attributes, respectively.

Accessing Domain Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Block factories may assume that all the publicly tunable parameters of
the domain have been set, such as the `mb->bc` boundary condition
parameters. These can be used to set the block mappings. As block
factories are run at the beginning of the `mrc_domain` setup process
one *cannot* assume private members of the domain have been
initialized.

