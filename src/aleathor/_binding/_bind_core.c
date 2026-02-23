// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Type lifecycle (new/init/dealloc), property getters,
 *           point queries, cell operations, universe operations.
 */

/* ============================================================================
 * Type Lifecycle
 * ============================================================================ */

static void AleaTHORSystem_dealloc(AleaTHORSystemObject* self) {
    if (self->sys && self->owns_sys) {
        alea_destroy(self->sys);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORSystem_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    AleaTHORSystemObject* self = (AleaTHORSystemObject*)type->tp_alloc(type, 0);
    if (self) {
        self->sys = NULL;
        self->owns_sys = 1;
    }
    return (PyObject*)self;
}

static int AleaTHORSystem_init(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "", kwlist)) {
        return -1;
    }

    self->sys = alea_create();
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create CSG system");
        return -1;
    }
    self->owns_sys = 1;
    return 0;
}

/* ============================================================================
 * Property Getters
 * ============================================================================ */

static PyObject* AleaTHORSystem_get_cell_count(AleaTHORSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_cell_count(self->sys));
}

static PyObject* AleaTHORSystem_get_surface_count(AleaTHORSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_surface_count(self->sys));
}

static PyObject* AleaTHORSystem_get_universe_count(AleaTHORSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_universe_count(self->sys));
}

static PyObject* AleaTHORSystem_get_spatial_index_instance_count(AleaTHORSystemObject* self, void* closure) {
    (void)closure;
    if (!self->sys) {
        return PyLong_FromLong(0);
    }
    return PyLong_FromSize_t(alea_spatial_index_instance_count(self->sys));
}

/* ============================================================================
 * Point Queries
 * ============================================================================ */

static PyObject* AleaTHORSystem_find_cell(AleaTHORSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int cell_idx = alea_find_cell(self->sys, x, y, z);
    if (cell_idx < 0) {
        Py_RETURN_NONE;  /* Point in void or outside */
    }

    int material_id = alea_material_at(self->sys, x, y, z);
    return Py_BuildValue("(ii)", cell_idx, material_id);
}

static PyObject* AleaTHORSystem_point_inside(AleaTHORSystemObject* self, PyObject* args) {
    unsigned long node_id;
    double x, y, z;

    if (!PyArg_ParseTuple(args, "kddd", &node_id, &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    bool inside = alea_point_inside(self->sys, (alea_node_id_t)node_id, x, y, z);
    return PyBool_FromLong(inside);
}

static PyObject* AleaTHORSystem_material_at(AleaTHORSystemObject* self, PyObject* args) {
    double x, y, z;

    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_material_id_t mat = alea_material_at(self->sys, x, y, z);
    if (mat == ALEA_MATERIAL_NONE) {
        Py_RETURN_NONE;
    }
    return PyLong_FromUnsignedLong(mat);
}

/* ============================================================================
 * Cell Operations
 * ============================================================================ */

static PyObject* AleaTHORSystem_get_cell(AleaTHORSystemObject* self, PyObject* args) {
    int cell_id;

    if (!PyArg_ParseTuple(args, "i", &cell_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_cell_info_t info;
    if (alea_cell_find_info(self->sys, cell_id, &info) < 0) {
        PyErr_Format(PyExc_KeyError, "Cell %d not found", cell_id);
        return NULL;
    }

    PyObject* dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:i, s:k, s:(dddddd)}",
        "cell_id", info.cell_id,
        "material_id", info.material_id,
        "density", info.density,
        "is_mass_density", PyBool_FromLong(info.is_mass_density),
        "universe_id", info.universe_id,
        "fill_universe", info.fill_universe,
        "fill_transform", info.fill_transform,
        "root_node", (unsigned long)info.root,
        "bbox", info.bbox.min_x, info.bbox.max_x,
                info.bbox.min_y, info.bbox.max_y,
                info.bbox.min_z, info.bbox.max_z);
    if (!dict) return NULL;

    if (add_lattice_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    return dict;
}

static PyObject* AleaTHORSystem_get_cell_by_index(AleaTHORSystemObject* self, PyObject* args) {
    size_t index;

    if (!PyArg_ParseTuple(args, "n", &index)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_cell_info_t info;
    if (alea_cell_get_info(self->sys, index, &info) < 0) {
        PyErr_Format(PyExc_IndexError, "Cell index %zu out of range", index);
        return NULL;
    }

    PyObject* dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:i, s:k, s:(dddddd)}",
        "cell_id", info.cell_id,
        "material_id", info.material_id,
        "density", info.density,
        "is_mass_density", PyBool_FromLong(info.is_mass_density),
        "universe_id", info.universe_id,
        "fill_universe", info.fill_universe,
        "fill_transform", info.fill_transform,
        "root_node", (unsigned long)info.root,
        "bbox", info.bbox.min_x, info.bbox.max_x,
                info.bbox.min_y, info.bbox.max_y,
                info.bbox.min_z, info.bbox.max_z);
    if (!dict) return NULL;

    if (add_lattice_fields(dict, &info) < 0) {
        Py_DECREF(dict);
        return NULL;
    }
    return dict;
}

static PyObject* AleaTHORSystem_get_cells(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_cell_count(self->sys);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        alea_cell_info_t info;
        if (alea_cell_get_info(self->sys, i, &info) < 0) {
            Py_DECREF(list);
            PyErr_SetString(PyExc_RuntimeError, "Failed to get cell info");
            return NULL;
        }

        PyObject* cell_dict = Py_BuildValue("{s:i, s:i, s:d, s:N, s:i, s:i, s:k}",
            "cell_id", info.cell_id,
            "material_id", info.material_id,
            "density", info.density,
            "is_mass_density", PyBool_FromLong(info.is_mass_density),
            "universe_id", info.universe_id,
            "fill_universe", info.fill_universe,
            "root_node", (unsigned long)info.root);

        if (!cell_dict) {
            Py_DECREF(list);
            return NULL;
        }

        if (add_lattice_fields(cell_dict, &info) < 0) {
            Py_DECREF(cell_dict);
            Py_DECREF(list);
            return NULL;
        }

        PyList_SET_ITEM(list, i, cell_dict);
    }

    return list;
}

/* ============================================================================
 * Universe Operations
 * ============================================================================ */

static PyObject* AleaTHORSystem_build_universe_index(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (alea_build_universe_index(self->sys) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_flatten_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id = 0;

    if (!PyArg_ParseTuple(args, "|i", &universe_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_flatten(self->sys, universe_id);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return PyLong_FromLong(result);
}

static PyObject* AleaTHORSystem_build_spatial_index(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_build_spatial_index(self->sys);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_get_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id;

    if (!PyArg_ParseTuple(args, "i", &universe_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int idx = alea_universe_find(self->sys, universe_id);
    if (idx < 0) {
        PyErr_Format(PyExc_KeyError, "Universe %d not found", universe_id);
        return NULL;
    }

    int uid;
    size_t cell_count;
    alea_bbox_t bbox;
    if (alea_universe_get(self->sys, (size_t)idx, &uid, &cell_count, &bbox) < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to get universe %d info", universe_id);
        return NULL;
    }

    return Py_BuildValue("{s:i, s:n, s:(dddddd)}",
        "universe_id", uid,
        "cell_count", (Py_ssize_t)cell_count,
        "bbox", bbox.min_x, bbox.max_x,
                bbox.min_y, bbox.max_y,
                bbox.min_z, bbox.max_z);
}
