======================
Parameter Interface
======================

.. c:type:: param_obj

   .. c:function:: void CLASS_set_param_obj(obj_type *obj, const char *name, void *val)

   .. c:function:: void *CLASS_get_param_obj(obj_type *obj, const char *name)


.. c:type:: param_int3

   .. c:function:: void CLASS_set_param_int3(obj_type *obj, const char *name, const int val[3])

   .. c:function:: void CLASS_get_param_int3(obj_type *obj, const char *name, int *pval)


.. c:type:: param_string

   .. c:function:: void CLASS_set_param_string(obj_type *obj, const char *name, const char *val)

   .. c:function:: void CLASS_get_param_string(obj_type *obj, const char *name, const char **val)
