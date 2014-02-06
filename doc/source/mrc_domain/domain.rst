====================
Domain
====================


Class Overview
================

In libmrc the simulation space is defined via an mrc_domain object.

.. _mrc_domain_types:

Types
------

.. toctree::
   :maxdepth: 2

   domain_mb



User Interface
==================

.. c:type:: struct mrc_domain

   The basic domain class doesn't have any parameters that need to be
   set or publically visible members which can be accessed by the
   standard `get_param` interface. Most generic interaction is handled
   through the helper functions below.

.. note:: While the main :c:type:`mrc_domain` class doesn't have
	  params or members, the individual types do. Consult the
	  individual :ref:`type documentation <mrc_domain_types>`
	  for the user interface.


Associated Functions
----------------------

.. c:function:: void mrc_domain_get_global_dims(struct mrc_domain *domain, int *dims)

   Get the global simulation size in each of the three spatial dimensions.

   :param int* dims: A pointer to store store the found dimensions.
 
.. warning:: A 3D cartesian volume is ill-defined for :c:type:`multi-block <mrc_domain_mb>` domains and will
	     trip an assert.


.. c:function:: void mrc_domain_get_bc(struct mrc_domain *domain, int *bc)

.. c:function:: void mrc_domain_get_nr_procs(struct mrc_domain *domain, int *nr_procs)

.. c:function:: void mrc_domain_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)

.. c:function:: void mrc_domain_get_global_patch_info(struct mrc_domain *domain, int gpatch, struct mrc_patch_info *info)

.. c:function:: void mrc_domain_get_local_patch_info(struct mrc_domain *domain, int patch, struct mrc_patch_info *info)

.. c:function:: void mrc_domain_get_level_idx3_patch_info(struct mrc_domain *domain, int level, int idx[3], struct mrc_patch_info *info)

.. c:function:: void mrc_domain_get_nr_levels(struct mrc_domain *domain, int *p_nr_levels)

.. c:function:: void mrc_domain_plot(struct mrc_domain *domain)

.. c:function:: int  mrc_domain_get_neighbor_rank(struct mrc_domain *domain, int shift[3])

.. c:function:: void mrc_domain_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3], int *nei_rank, int *nei_patch)

.. c:function:: void mrc_domain_add_patch(struct mrc_domain *domain, int l, int idx3[3])

.. c:function:: struct mrc_patch* mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches)

.. c:function:: struct mrc_crds* mrc_domain_get_crds(struct mrc_domain *domain)

.. c:function:: struct mrc_ddc* mrc_domain_get_ddc(struct mrc_domain *domain)

.. c:function:: struct mrc_fld* mrc_domain_f1_create(struct mrc_domain *domain)

.. c:function:: struct mrc_m3* mrc_domain_m3_create(struct mrc_domain *domain)

.. c:function:: struct mrc_fld* mrc_domain_m1_create(struct mrc_domain *domain)

.. c:function:: struct mrc_fld* mrc_domain_fld_create(struct mrc_domain *domain, int sw, const char *comps)

.. c:function:: struct mrc_ddc* mrc_domain_create_ddc(struct mrc_domain *domain)

Writing A Subclass
===================

Don't do it. Seriously. You're asking for a world of trouble. I
suppose it's possible though. I mean, I sort of did it. So I guess we
should probably write something here. I just don't want to be the guy
sitting in the bed of a pickup encouraging my buddy to piss on the
electrified fence.
