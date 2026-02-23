// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Primitive creation, boolean operations, surface creation,
 *           cell registration (add_cell).
 */

/* ============================================================================
 * AleaTHORSystem Methods - Primitive Creation
 * ============================================================================ */

static PyObject* AleaTHORSystem_create_plane(AleaTHORSystemObject* self, PyObject* args) {
    double a, b, c, d;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddi", &a, &b, &c, &d, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_plane_surface(self->sys, 0, a, b, c, d);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_sphere(AleaTHORSystemObject* self, PyObject* args) {
    double cx, cy, cz, radius;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddi", &cx, &cy, &cz, &radius, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sphere_surface(self->sys, 0, cx, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_box(AleaTHORSystemObject* self, PyObject* args) {
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int sense;
    if (!PyArg_ParseTuple(args, "ddddddi", &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_surface(self->sys, 0, xmin, xmax, ymin, ymax, zmin, zmax);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_cylinder_z(AleaTHORSystemObject* self, PyObject* args) {
    double cx, cy, radius;
    int sense;
    if (!PyArg_ParseTuple(args, "dddi", &cx, &cy, &radius, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_z_surface(self->sys, 0, cx, cy, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder");
        return NULL;
    }
    alea_node_id_t node = alea_halfspace(self->sys, idx, sense);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder half-space");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

/* ============================================================================
 * AleaTHORSystem Methods - Boolean Operations
 * ============================================================================ */

static PyObject* AleaTHORSystem_create_union(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_union(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create union");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_intersection(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_intersection(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create intersection");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_difference(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long a, b;
    if (!PyArg_ParseTuple(args, "kk", &a, &b)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_difference(self->sys, (alea_node_id_t)a, (alea_node_id_t)b);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create difference");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_complement(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long a;
    if (!PyArg_ParseTuple(args, "k", &a)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t node = alea_complement(self->sys, (alea_node_id_t)a);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create complement");
        return NULL;
    }
    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORSystem_create_union_many(AleaTHORSystemObject* self, PyObject* args) {
    PyObject* nodes_list;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &nodes_list)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    Py_ssize_t count = PyList_Size(nodes_list);
    if (count < 1) {
        PyErr_SetString(PyExc_ValueError, "Need at least one node");
        return NULL;
    }

    alea_node_id_t* nodes = malloc(count * sizeof(alea_node_id_t));
    if (!nodes) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = PyList_GetItem(nodes_list, i);
        if (!PyLong_Check(item)) {
            free(nodes);
            PyErr_SetString(PyExc_TypeError, "All items must be integers (node IDs)");
            return NULL;
        }
        nodes[i] = (alea_node_id_t)PyLong_AsUnsignedLong(item);
    }

    alea_node_id_t result = alea_union_n(self->sys, nodes, count);
    free(nodes);

    if (result == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create union");
        return NULL;
    }
    return PyLong_FromUnsignedLong(result);
}

static PyObject* AleaTHORSystem_create_intersection_many(AleaTHORSystemObject* self, PyObject* args) {
    PyObject* nodes_list;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &nodes_list)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    Py_ssize_t count = PyList_Size(nodes_list);
    if (count < 1) {
        PyErr_SetString(PyExc_ValueError, "Need at least one node");
        return NULL;
    }

    alea_node_id_t* nodes = malloc(count * sizeof(alea_node_id_t));
    if (!nodes) {
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = PyList_GetItem(nodes_list, i);
        if (!PyLong_Check(item)) {
            free(nodes);
            PyErr_SetString(PyExc_TypeError, "All items must be integers (node IDs)");
            return NULL;
        }
        nodes[i] = (alea_node_id_t)PyLong_AsUnsignedLong(item);
    }

    alea_node_id_t result = alea_intersection_n(self->sys, nodes, count);
    free(nodes);

    if (result == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create intersection");
        return NULL;
    }
    return PyLong_FromUnsignedLong(result);
}

/* ============================================================================
 * AleaTHORSystem Methods - Surface Creation (with automatic registration)
 *
 * These functions use alea_*_surface() which creates surfaces with both
 * positive and negative halfspace nodes, properly registered for raycast.
 * ============================================================================ */

static PyObject* AleaTHORSystem_sphere_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, radius;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sphere_surface(self->sys, surface_id, cx, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create sphere surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cylinder_z_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cx, &cy, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_z_surface(self->sys, surface_id, cx, cy, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_box_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    if (!PyArg_ParseTuple(args, "idddddd", &surface_id, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_surface(self->sys, surface_id, xmin, xmax, ymin, ymax, zmin, zmax);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create box surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_plane_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double a, b, c, d;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &a, &b, &c, &d)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_plane_surface(self->sys, surface_id, a, b, c, d);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create plane surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cylinder_x_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cy, cz, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cy, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_x_surface(self->sys, surface_id, cy, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cylinder_y_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cz, radius;
    if (!PyArg_ParseTuple(args, "iddd", &surface_id, &cx, &cz, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cylinder_y_surface(self->sys, surface_id, cx, cz, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cylinder surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cone_z_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_z_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cone_x_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_x_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_cone_y_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, t_squared;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &t_squared)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_cone_y_surface(self->sys, surface_id, cx, cy, cz, t_squared);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create cone surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_torus_z_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_z_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_torus_x_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_x_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_torus_y_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, major_radius, minor_radius;
    if (!PyArg_ParseTuple(args, "iddddd", &surface_id, &cx, &cy, &cz, &major_radius, &minor_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_torus_y_surface(self->sys, surface_id, cx, cy, cz, major_radius, minor_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create torus surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_quadric_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double A, B, C, D, E, F, G, H, I, J;
    if (!PyArg_ParseTuple(args, "idddddddddd", &surface_id, &A, &B, &C, &D, &E, &F, &G, &H, &I, &J)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_quadric_surface(self->sys, surface_id, A, B, C, D, E, F, G, H, I, J);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create quadric surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_rcc_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z, radius;
    if (!PyArg_ParseTuple(args, "iddddddd", &surface_id, &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z, &radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rcc_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z, radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create RCC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_box_general_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double corner_x, corner_y, corner_z;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &corner_x, &corner_y, &corner_z,
                          &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z,
                          &v3_x, &v3_y, &v3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_box_general_surface(self->sys, surface_id, corner_x, corner_y, corner_z,
                                           v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create general box surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_sph_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double cx, cy, cz, r;
    if (!PyArg_ParseTuple(args, "idddd", &surface_id, &cx, &cy, &cz, &r)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_sph_surface(self->sys, surface_id, cx, cy, cz, r);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create SPH surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_trc_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z, base_radius, top_radius;
    if (!PyArg_ParseTuple(args, "idddddddd", &surface_id, &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z, &base_radius, &top_radius)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_trc_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z, base_radius, top_radius);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create TRC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_ell_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, major_axis_len;
    if (!PyArg_ParseTuple(args, "iddddddd", &surface_id, &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z, &major_axis_len)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_ell_surface(self->sys, surface_id, v1_x, v1_y, v1_z,
                                   v2_x, v2_y, v2_z, major_axis_len);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create ELL surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_rec_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z;
    double axis1_x, axis1_y, axis1_z, axis2_x, axis2_y, axis2_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z,
                          &axis1_x, &axis1_y, &axis1_z,
                          &axis2_x, &axis2_y, &axis2_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rec_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z,
                                   axis1_x, axis1_y, axis1_z,
                                   axis2_x, axis2_y, axis2_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create REC surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_wed_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double vertex_x, vertex_y, vertex_z;
    double v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddd", &surface_id,
                          &vertex_x, &vertex_y, &vertex_z,
                          &v1_x, &v1_y, &v1_z,
                          &v2_x, &v2_y, &v2_z,
                          &v3_x, &v3_y, &v3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_wed_surface(self->sys, surface_id, vertex_x, vertex_y, vertex_z,
                                   v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create WED surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_rhp_surface(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    double base_x, base_y, base_z, height_x, height_y, height_z;
    double r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, r3_x, r3_y, r3_z;
    if (!PyArg_ParseTuple(args, "idddddddddddddddd", &surface_id,
                          &base_x, &base_y, &base_z,
                          &height_x, &height_y, &height_z,
                          &r1_x, &r1_y, &r1_z,
                          &r2_x, &r2_y, &r2_z,
                          &r3_x, &r3_y, &r3_z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_rhp_surface(self->sys, surface_id, base_x, base_y, base_z,
                                   height_x, height_y, height_z,
                                   r1_x, r1_y, r1_z,
                                   r2_x, r2_y, r2_z,
                                   r3_x, r3_y, r3_z);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create RHP surface");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(ikk)", idx, (unsigned long)pos_node, (unsigned long)neg_node);
}

static PyObject* AleaTHORSystem_get_surface_nodes(AleaTHORSystemObject* self, PyObject* args) {
    int surface_idx;
    if (!PyArg_ParseTuple(args, "i", &surface_idx)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (surface_idx < 0 || (size_t)surface_idx >= alea_surface_count(self->sys)) {
        PyErr_SetString(PyExc_IndexError, "Surface index out of range");
        return NULL;
    }

    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, surface_idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return Py_BuildValue("(kk)", (unsigned long)pos_node, (unsigned long)neg_node);
}

/* ============================================================================
 * AleaTHORSystem Methods - Cell Registration
 * ============================================================================ */

static PyObject* AleaTHORSystem_add_cell(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    int cell_id;
    unsigned long root_node;
    int material_id = 0;
    double density = 0.0;
    int universe_id = 0;
    static char* kwlist[] = {"cell_id", "root_node", "material_id", "density", "universe_id", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ik|idi", kwlist,
                                     &cell_id, &root_node, &material_id, &density, &universe_id)) {
        return NULL;
    }
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_add_cell(self->sys, cell_id, (alea_node_id_t)root_node, material_id, density, universe_id);
    if (idx < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to add cell");
        return NULL;
    }
    return PyLong_FromLong(idx);
}
