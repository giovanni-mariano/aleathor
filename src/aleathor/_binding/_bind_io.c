// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Export (MCNP, OpenMC), merge, extract/filter operations,
 *           material operations, renumbering, split/expand.
 */

/* ============================================================================
 * AleaTHORSystem Methods - Export
 * ============================================================================ */

static PyObject* AleaTHORSystem_export_mcnp(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    int deduplicate = 1;
    int universe_depth = -1;  /* -1 = all universes */
    int fill_depth = 0;       /* 0 = no expansion */
    static char* kwlist[] = {"filename", "deduplicate", "universe_depth", "fill_depth", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|pii", kwlist,
            &filename, &deduplicate, &universe_depth, &fill_depth)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Save original config and apply export settings */
    alea_config_t orig = alea_get_config(self->sys);
    alea_config_t cfg = orig;
    cfg.dedup = deduplicate;
    cfg.universe_depth = universe_depth;
    cfg.fill_depth = fill_depth;
    alea_set_config(self->sys, &cfg);

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_export_mcnp(self->sys, filename);
    Py_END_ALLOW_THREADS

    /* Restore original config */
    alea_set_config(self->sys, &orig);

    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_Format(PyExc_IOError, "Failed to export to %s: %s", filename, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_export_openmc(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    static char* kwlist[] = {"filename", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &filename)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_export_openmc(self->sys, filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_Format(PyExc_IOError, "Failed to export to %s: %s", filename, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

/* ============================================================================
 * AleaTHORSystem Methods - Merge
 * ============================================================================ */

static PyObject* AleaTHORSystem_merge(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    AleaTHORSystemObject* other;
    int id_offset = 0;
    static char* kwlist[] = {"other", "id_offset", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|i", kwlist, &AleaTHORSystemType, &other, &id_offset)) {
        return NULL;
    }

    if (!self->sys || !other->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (alea_merge(self->sys, other->sys, id_offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

/* ============================================================================
 * AleaTHORSystem Methods - Extract / Filter
 * ============================================================================ */

static PyObject* AleaTHORSystem_get_cells_by_material(AleaTHORSystemObject* self, PyObject* args) {
    int material_id;
    if (!PyArg_ParseTuple(args, "i", &material_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* First call to get count */
    size_t count = alea_get_cells_by_material(self->sys, material_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_by_material(self->sys, material_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* AleaTHORSystem_get_cells_by_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_get_cells_by_universe(self->sys, universe_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_by_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* AleaTHORSystem_get_cells_filling_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    size_t count = alea_get_cells_filling_universe(self->sys, universe_id, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_filling_universe(self->sys, universe_id, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* AleaTHORSystem_get_cells_in_bbox(AleaTHORSystemObject* self, PyObject* args) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    if (!PyArg_ParseTuple(args, "dddddd", &x_min, &x_max, &y_min, &y_max, &z_min, &z_max)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_bbox_t bbox = { x_min, x_max, y_min, y_max, z_min, z_max };
    size_t count = alea_get_cells_in_bbox(self->sys, &bbox, NULL, 0);

    if (count == 0) {
        return PyList_New(0);
    }

    int* indices = malloc(count * sizeof(int));
    if (!indices) return PyErr_NoMemory();

    alea_get_cells_in_bbox(self->sys, &bbox, indices, count);

    PyObject* result = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(indices[i]));
    }

    free(indices);
    return result;
}

static PyObject* AleaTHORSystem_extract_universe(AleaTHORSystemObject* self, PyObject* args) {
    int universe_id;
    if (!PyArg_ParseTuple(args, "i", &universe_id)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_system_t* new_sys = alea_extract_universe(self->sys, universe_id);
    if (!new_sys) {
        PyErr_Format(PyExc_RuntimeError, "Failed to extract universe: %s", alea_error());
        return NULL;
    }

    /* Create new Python object wrapping the extracted system */
    AleaTHORSystemObject* new_obj = (AleaTHORSystemObject*)AleaTHORSystem_new(&AleaTHORSystemType, NULL, NULL);
    if (!new_obj) {
        alea_destroy(new_sys);
        return NULL;
    }

    new_obj->sys = new_sys;
    new_obj->owns_sys = 1;

    return (PyObject*)new_obj;
}

static PyObject* AleaTHORSystem_extract_region(AleaTHORSystemObject* self, PyObject* args) {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    if (!PyArg_ParseTuple(args, "dddddd", &x_min, &x_max, &y_min, &y_max, &z_min, &z_max)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_bbox_t bbox = { x_min, x_max, y_min, y_max, z_min, z_max };
    alea_system_t* new_sys = alea_extract_region(self->sys, &bbox);
    if (!new_sys) {
        PyErr_Format(PyExc_RuntimeError, "Failed to extract region: %s", alea_error());
        return NULL;
    }

    AleaTHORSystemObject* new_obj = (AleaTHORSystemObject*)AleaTHORSystem_new(&AleaTHORSystemType, NULL, NULL);
    if (!new_obj) {
        alea_destroy(new_sys);
        return NULL;
    }

    new_obj->sys = new_sys;
    new_obj->owns_sys = 1;

    return (PyObject*)new_obj;
}

/* ============================================================================
 * AleaTHORSystem Methods - Material Operations
 * ============================================================================ */

static PyObject* AleaTHORSystem_create_mixture(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* mat_list;
    PyObject* frac_list;
    int new_mat_id = 0;
    static char* kwlist[] = {"material_ids", "fractions", "new_id", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|i", kwlist, &mat_list, &frac_list, &new_mat_id)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (!PyList_Check(mat_list) || !PyList_Check(frac_list)) {
        PyErr_SetString(PyExc_TypeError, "material_ids and fractions must be lists");
        return NULL;
    }

    Py_ssize_t count = PyList_Size(mat_list);
    if (count != PyList_Size(frac_list)) {
        PyErr_SetString(PyExc_ValueError, "material_ids and fractions must have same length");
        return NULL;
    }

    if (count == 0) {
        PyErr_SetString(PyExc_ValueError, "At least one material required");
        return NULL;
    }

    int* mat_ids = malloc(count * sizeof(int));
    double* fractions = malloc(count * sizeof(double));
    if (!mat_ids || !fractions) {
        free(mat_ids);
        free(fractions);
        return PyErr_NoMemory();
    }

    for (Py_ssize_t i = 0; i < count; i++) {
        mat_ids[i] = (int)PyLong_AsLong(PyList_GetItem(mat_list, i));
        fractions[i] = PyFloat_AsDouble(PyList_GetItem(frac_list, i));
    }

    if (PyErr_Occurred()) {
        free(mat_ids);
        free(fractions);
        return NULL;
    }

    int result = alea_create_mixture(self->sys, mat_ids, fractions, (size_t)count, new_mat_id);

    free(mat_ids);
    free(fractions);

    if (result < 0) {
        PyErr_Format(PyExc_RuntimeError, "Failed to create mixture: %s", alea_error());
        return NULL;
    }

    return PyLong_FromLong(result);
}

/* ============================================================================
 * AleaTHORSystem Methods - Renumbering
 * ============================================================================ */

static PyObject* AleaTHORSystem_renumber_cells(AleaTHORSystemObject* self, PyObject* args) {
    int start_id;
    if (!PyArg_ParseTuple(args, "i", &start_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_renumber_cells(self->sys, start_id);
    if (result < 0) { PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL; }
    return PyLong_FromLong(result);
}

static PyObject* AleaTHORSystem_renumber_surfaces(AleaTHORSystemObject* self, PyObject* args) {
    int start_id;
    if (!PyArg_ParseTuple(args, "i", &start_id)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_renumber_surfaces(self->sys, start_id);
    if (result < 0) { PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL; }
    return PyLong_FromLong(result);
}

static PyObject* AleaTHORSystem_offset_cell_ids(AleaTHORSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_cell_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_offset_surface_ids(AleaTHORSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_surface_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_offset_material_ids(AleaTHORSystemObject* self, PyObject* args) {
    int offset;
    if (!PyArg_ParseTuple(args, "i", &offset)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    if (alea_offset_material_ids(self->sys, offset) < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error()); return NULL;
    }
    Py_RETURN_NONE;
}

/* ============================================================================
 * AleaTHORSystem Methods - Split / Expand
 * ============================================================================ */

static PyObject* AleaTHORSystem_split_union_cells(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_split_union_cells(self->sys);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

static PyObject* AleaTHORSystem_expand_macrobodies(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    int result = alea_expand_macrobodies_in_system(self->sys);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}
