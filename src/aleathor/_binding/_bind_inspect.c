// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Cell lookup, find all cells at point, set fill, cell/surface/node info,
 *           node inspection, bbox tightening, CSG simplification, numerical bbox,
 *           primitive data, CSG node tree.
 */

/* ============================================================================
 * AleaTHORSystem Methods - Cell Lookup by ID
 * ============================================================================ */

static PyObject* AleaTHORSystem_cell_find(AleaTHORSystemObject* self, PyObject* args) {
    int cell_id;

    if (!PyArg_ParseTuple(args, "i", &cell_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int idx = alea_cell_find(self->sys, cell_id);
    if (idx < 0) {
        Py_RETURN_NONE;
    }
    return PyLong_FromLong(idx);
}

/* ============================================================================
 * AleaTHORSystem Methods - Find All Cells at Point (Hierarchy)
 * ============================================================================ */

static PyObject* AleaTHORSystem_find_all_cells(AleaTHORSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Allocate hits on stack - 64 levels should be more than enough */
    alea_cell_hit_t hits[64];
    int nhits = alea_find_all_cells(self->sys, x, y, z, hits, 64);
    if (nhits < 0) {
        PyErr_SetString(PyExc_RuntimeError, "find_all_cells failed");
        return NULL;
    }

    PyObject* list = PyList_New(nhits);
    if (!list) return NULL;

    for (int i = 0; i < nhits; i++) {
        PyObject* hit = Py_BuildValue("{s:i, s:i, s:i, s:i, s:i, s:i, s:d, s:d, s:d}",
            "cell_id", hits[i].cell_id,
            "cell_index", hits[i].cell_index,
            "material_id", hits[i].material_id,
            "universe_id", hits[i].universe_id,
            "fill_universe", hits[i].fill_universe,
            "depth", hits[i].depth,
            "local_x", hits[i].local_x,
            "local_y", hits[i].local_y,
            "local_z", hits[i].local_z);
        if (!hit) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, hit);
    }

    return list;
}

/* ============================================================================
 * AleaTHORSystem Methods - Set Fill
 * ============================================================================ */

static PyObject* AleaTHORSystem_set_fill(AleaTHORSystemObject* self, PyObject* args) {
    int cell_index, fill_universe, transform = 0;
    if (!PyArg_ParseTuple(args, "ii|i", &cell_index, &fill_universe, &transform)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_set_fill(self->sys, cell_index, fill_universe, transform) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to set fill: %s", alea_error());
        return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * AleaTHORSystem Methods - Cell/Surface/Node Info
 * ============================================================================ */

static PyObject* AleaTHORSystem_get_cell_id(AleaTHORSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int cell_id = alea_get_cell_id(self->sys, cell_index);
    if (cell_id < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %d out of range", cell_index);
        return NULL;
    }
    return PyLong_FromLong(cell_id);
}

static PyObject* AleaTHORSystem_surface_find(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id;
    if (!PyArg_ParseTuple(args, "i", &surface_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_surface_find(self->sys, surface_id);
    if (idx < 0) Py_RETURN_NONE;
    return PyLong_FromLong(idx);
}

static PyObject* AleaTHORSystem_surface_node(AleaTHORSystemObject* self, PyObject* args) {
    int surface_id, sense;
    if (!PyArg_ParseTuple(args, "ii", &surface_id, &sense)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int idx = alea_surface_find(self->sys, surface_id);
    if (idx < 0) Py_RETURN_NONE;
    alea_node_id_t pos_node, neg_node;
    alea_surface_get(self->sys, idx, NULL, NULL, &pos_node, &neg_node, NULL);
    return PyLong_FromUnsignedLong(sense > 0 ? pos_node : neg_node);
}

static PyObject* AleaTHORSystem_cells_in_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    /* First call to get count */
    int count = alea_cells_in_universe(self->sys, universe_id, NULL, 0);
    if (count <= 0) return PyList_New(0);

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();
    alea_cells_in_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }
    free(indices);
    return result;
}

static PyObject* AleaTHORSystem_find_cell_at(AleaTHORSystemObject* self, PyObject* args) {
    double x, y, z;
    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int cell_id, material;
    int result = alea_find_cell_at(self->sys, x, y, z, &cell_id, &material);
    if (result < 0) Py_RETURN_NONE;
    return Py_BuildValue("(ii)", cell_id, material);
}

/* ============================================================================
 * AleaTHORSystem Methods - Node Inspection
 * ============================================================================ */

static PyObject* AleaTHORSystem_node_primitive_type(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_primitive_type_t ptype = alea_node_primitive_type(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong((int)ptype);
}

static PyObject* AleaTHORSystem_node_primitive_id(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_primitive_id_t pid = alea_node_primitive_id(self->sys, (alea_node_id_t)node_id);
    if (pid == ALEA_PRIMITIVE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(pid);
}

static PyObject* AleaTHORSystem_node_sense(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int sense = alea_node_sense(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong(sense);
}

static PyObject* AleaTHORSystem_node_surface_id(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int sid = alea_node_surface_id(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong(sid);
}

static PyObject* AleaTHORSystem_node_operation(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_operation_t op = alea_node_operation(self->sys, (alea_node_id_t)node_id);
    return PyLong_FromLong((int)op);
}

static PyObject* AleaTHORSystem_node_left(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t left = alea_node_left(self->sys, (alea_node_id_t)node_id);
    if (left == ALEA_NODE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(left);
}

static PyObject* AleaTHORSystem_node_right(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_node_id_t right = alea_node_right(self->sys, (alea_node_id_t)node_id);
    if (right == ALEA_NODE_ID_INVALID) Py_RETURN_NONE;
    return PyLong_FromUnsignedLong(right);
}

/* ============================================================================
 * AleaTHORSystem Methods - BBox Tightening
 * ============================================================================ */

static PyObject* AleaTHORSystem_tighten_cell_bbox(AleaTHORSystemObject* self, PyObject* args) {
    size_t cell_index;
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "n|d", &cell_index, &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_bbox_t out;
    if (alea_tighten_cell_bbox(self->sys, cell_index, tol, &out) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return Py_BuildValue("(dddddd)", out.min_x, out.max_x, out.min_y, out.max_y, out.min_z, out.max_z);
}

static PyObject* AleaTHORSystem_tighten_all_bboxes(AleaTHORSystemObject* self, PyObject* args) {
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "|d", &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_tighten_all_bboxes(self->sys, tol);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

/* ============================================================================
 * AleaTHORSystem Methods - CSG Simplification
 * ============================================================================ */

static PyObject* AleaTHORSystem_flatten_all_cells(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_simplify_stats_t stats;
    memset(&stats, 0, sizeof(stats));

    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    alea_flatten_all_cells(self->sys, &stats);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    return Py_BuildValue(
        "{s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n, s:n}",
        "nodes_before",             (Py_ssize_t)stats.nodes_before,
        "nodes_after",              (Py_ssize_t)stats.nodes_after,
        "complements_eliminated",   (Py_ssize_t)stats.complements_eliminated,
        "double_negations",         (Py_ssize_t)stats.double_negations,
        "idempotent_reductions",    (Py_ssize_t)stats.idempotent_reductions,
        "absorption_reductions",    (Py_ssize_t)stats.absorption_reductions,
        "subtrees_deduplicated",    (Py_ssize_t)stats.subtrees_deduplicated,
        "cell_complements_expanded",(Py_ssize_t)stats.cell_complements_expanded,
        "contradictions_found",     (Py_ssize_t)stats.contradictions_found,
        "tautologies_found",        (Py_ssize_t)stats.tautologies_found,
        "empty_cells_removed",      (Py_ssize_t)stats.empty_cells_removed,
        "union_branches_absorbed",  (Py_ssize_t)stats.union_branches_absorbed,
        "union_common_factors",     (Py_ssize_t)stats.union_common_factors,
        "union_branches_subsumed",  (Py_ssize_t)stats.union_branches_subsumed
    );
}

/* ============================================================================
 * AleaTHORSystem Methods - Numerical BBox Tightening
 * ============================================================================ */

static PyObject* AleaTHORSystem_tighten_cell_bbox_numerical(AleaTHORSystemObject* self, PyObject* args) {
    int cell_index;
    if (!PyArg_ParseTuple(args, "i", &cell_index)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_tighten_cell_bbox_numerical(self->sys, cell_index);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

/* ============================================================================
 * AleaTHORSystem Methods - Primitive Data
 * ============================================================================ */

static PyObject* AleaTHORSystem_node_primitive_data(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_primitive_data_t data;
    int rc = alea_node_primitive_data(self->sys, (alea_node_id_t)node_id, &data);
    if (rc < 0) {
        PyErr_SetString(PyExc_ValueError, "Node is not a primitive or invalid");
        return NULL;
    }

    alea_primitive_type_t ptype = alea_node_primitive_type(self->sys, (alea_node_id_t)node_id);

    switch (ptype) {
    case ALEA_PRIMITIVE_PLANE:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "a", data.plane.a, "b", data.plane.b,
            "c", data.plane.c, "d", data.plane.d);

    case ALEA_PRIMITIVE_SPHERE:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.sphere.center_x,
            "center_y", data.sphere.center_y,
            "center_z", data.sphere.center_z,
            "radius", data.sphere.radius);

    case ALEA_PRIMITIVE_CYLINDER_X:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_y", data.cyl_x.center_y,
            "center_z", data.cyl_x.center_z,
            "radius", data.cyl_x.radius);

    case ALEA_PRIMITIVE_CYLINDER_Y:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.cyl_y.center_x,
            "center_z", data.cyl_y.center_z,
            "radius", data.cyl_y.radius);

    case ALEA_PRIMITIVE_CYLINDER_Z:
        return Py_BuildValue("{s:i, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.cyl_z.center_x,
            "center_y", data.cyl_z.center_y,
            "radius", data.cyl_z.radius);

    case ALEA_PRIMITIVE_CONE_X:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_x.apex_x,
            "apex_y", data.cone_x.apex_y,
            "apex_z", data.cone_x.apex_z,
            "tan_angle_sq", data.cone_x.tan_angle_sq,
            "sheet_selection", data.cone_x.sheet_selection);

    case ALEA_PRIMITIVE_CONE_Y:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_y.apex_x,
            "apex_y", data.cone_y.apex_y,
            "apex_z", data.cone_y.apex_z,
            "tan_angle_sq", data.cone_y.tan_angle_sq,
            "sheet_selection", data.cone_y.sheet_selection);

    case ALEA_PRIMITIVE_CONE_Z:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:i}",
            "type", (int)ptype,
            "apex_x", data.cone_z.apex_x,
            "apex_y", data.cone_z.apex_y,
            "apex_z", data.cone_z.apex_z,
            "tan_angle_sq", data.cone_z.tan_angle_sq,
            "sheet_selection", data.cone_z.sheet_selection);

    case ALEA_PRIMITIVE_RPP:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "min_x", data.box.min_x, "max_x", data.box.max_x,
            "min_y", data.box.min_y, "max_y", data.box.max_y,
            "min_z", data.box.min_z, "max_z", data.box.max_z);

    case ALEA_PRIMITIVE_QUADRIC: {
        PyObject* coeffs = PyList_New(10);
        if (!coeffs) return NULL;
        for (int i = 0; i < 10; i++) {
            PyList_SET_ITEM(coeffs, i, PyFloat_FromDouble(data.quadric.coeffs[i]));
        }
        PyObject* result = Py_BuildValue("{s:i, s:N}",
            "type", (int)ptype, "coeffs", coeffs);
        return result;
    }

    case ALEA_PRIMITIVE_TORUS_X:
    case ALEA_PRIMITIVE_TORUS_Y:
    case ALEA_PRIMITIVE_TORUS_Z:
        return Py_BuildValue("{s:i, s:i, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "axis", (int)data.torus.axis,
            "center_x", data.torus.center_x,
            "center_y", data.torus.center_y,
            "center_z", data.torus.center_z,
            "major_radius", data.torus.major_radius,
            "minor_radius", data.torus.minor_radius,
            "axial_semiwidth_B", data.torus.axial_semiwidth_B);

    case ALEA_PRIMITIVE_RCC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rcc.base_x, "base_y", data.rcc.base_y, "base_z", data.rcc.base_z,
            "height_x", data.rcc.height_x, "height_y", data.rcc.height_y, "height_z", data.rcc.height_z,
            "radius", data.rcc.radius);

    case ALEA_PRIMITIVE_BOX:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "corner_x", data.box_general.corner_x,
            "corner_y", data.box_general.corner_y,
            "corner_z", data.box_general.corner_z,
            "v1_x", data.box_general.v1_x, "v1_y", data.box_general.v1_y, "v1_z", data.box_general.v1_z,
            "v2_x", data.box_general.v2_x, "v2_y", data.box_general.v2_y, "v2_z", data.box_general.v2_z,
            "v3_x", data.box_general.v3_x, "v3_y", data.box_general.v3_y, "v3_z", data.box_general.v3_z);

    case ALEA_PRIMITIVE_SPH:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "center_x", data.sph.center_x,
            "center_y", data.sph.center_y,
            "center_z", data.sph.center_z,
            "radius", data.sph.radius);

    case ALEA_PRIMITIVE_TRC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.trc.base_x, "base_y", data.trc.base_y, "base_z", data.trc.base_z,
            "height_x", data.trc.height_x, "height_y", data.trc.height_y, "height_z", data.trc.height_z,
            "base_radius", data.trc.base_radius, "top_radius", data.trc.top_radius);

    case ALEA_PRIMITIVE_ELL:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "v1_x", data.ell.v1_x, "v1_y", data.ell.v1_y, "v1_z", data.ell.v1_z,
            "v2_x", data.ell.v2_x, "v2_y", data.ell.v2_y, "v2_z", data.ell.v2_z,
            "major_axis_len", data.ell.major_axis_len);

    case ALEA_PRIMITIVE_REC:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rec.base_x, "base_y", data.rec.base_y, "base_z", data.rec.base_z,
            "height_x", data.rec.height_x, "height_y", data.rec.height_y, "height_z", data.rec.height_z,
            "axis1_x", data.rec.axis1_x, "axis1_y", data.rec.axis1_y, "axis1_z", data.rec.axis1_z,
            "axis2_x", data.rec.axis2_x, "axis2_y", data.rec.axis2_y, "axis2_z", data.rec.axis2_z);

    case ALEA_PRIMITIVE_WED:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "vertex_x", data.wed.vertex_x, "vertex_y", data.wed.vertex_y, "vertex_z", data.wed.vertex_z,
            "v1_x", data.wed.v1_x, "v1_y", data.wed.v1_y, "v1_z", data.wed.v1_z,
            "v2_x", data.wed.v2_x, "v2_y", data.wed.v2_y, "v2_z", data.wed.v2_z,
            "v3_x", data.wed.v3_x, "v3_y", data.wed.v3_y, "v3_z", data.wed.v3_z);

    case ALEA_PRIMITIVE_RHP:
        return Py_BuildValue("{s:i, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d, s:d}",
            "type", (int)ptype,
            "base_x", data.rhp.base_x, "base_y", data.rhp.base_y, "base_z", data.rhp.base_z,
            "height_x", data.rhp.height_x, "height_y", data.rhp.height_y, "height_z", data.rhp.height_z,
            "r1_x", data.rhp.r1_x, "r1_y", data.rhp.r1_y, "r1_z", data.rhp.r1_z,
            "r2_x", data.rhp.r2_x, "r2_y", data.rhp.r2_y, "r2_z", data.rhp.r2_z,
            "r3_x", data.rhp.r3_x, "r3_y", data.rhp.r3_y, "r3_z", data.rhp.r3_z);

    case ALEA_PRIMITIVE_ARB: {
        PyObject* corners = PyList_New(data.arb.num_corners);
        if (!corners) return NULL;
        for (int i = 0; i < data.arb.num_corners; i++) {
            PyObject* pt = Py_BuildValue("(ddd)",
                data.arb.corners[i][0], data.arb.corners[i][1], data.arb.corners[i][2]);
            if (!pt) { Py_DECREF(corners); return NULL; }
            PyList_SET_ITEM(corners, i, pt);
        }
        PyObject* faces = PyList_New(data.arb.num_faces);
        if (!faces) { Py_DECREF(corners); return NULL; }
        for (int i = 0; i < data.arb.num_faces; i++) {
            PyObject* face = Py_BuildValue("(iiii)",
                data.arb.faces[i][0], data.arb.faces[i][1],
                data.arb.faces[i][2], data.arb.faces[i][3]);
            if (!face) { Py_DECREF(corners); Py_DECREF(faces); return NULL; }
            PyList_SET_ITEM(faces, i, face);
        }
        PyObject* result = Py_BuildValue("{s:i, s:N, s:N, s:i, s:i}",
            "type", (int)ptype,
            "corners", corners, "faces", faces,
            "num_corners", data.arb.num_corners,
            "num_faces", data.arb.num_faces);
        return result;
    }

    default:
        return Py_BuildValue("{s:i}", "type", (int)ptype);
    }
}

static PyObject* AleaTHORSystem_node_tree(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    if (!PyArg_ParseTuple(args, "k", &node_id)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return build_node_tree(self->sys, (alea_node_id_t)node_id);
}
