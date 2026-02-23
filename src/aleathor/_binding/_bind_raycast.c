// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Raycast, cell-aware raycast, find overlaps,
 *           volume estimation, bounding sphere.
 */

/* ============================================================================
 * AleaTHORSystem Methods - Raycast
 * ============================================================================ */

static PyObject* AleaTHORSystem_raycast(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz;
    double t_max = 0.0;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|d", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) {
        return NULL;
    }
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_raycast_result_t* result = alea_raycast_result_create();
    if (!result) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate raycast result");
        return NULL;
    }

    int raycast_result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    raycast_result = alea_raycast(self->sys, ox, oy, oz, dx, dy, dz, t_max, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_raycast_result_free(result);
        return NULL;
    }

    if (raycast_result < 0) {
        alea_raycast_result_free(result);
        PyErr_SetString(PyExc_RuntimeError, "Raycast failed");
        return NULL;
    }

    size_t count = alea_raycast_segment_count(result);
    PyObject* segments = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id;
        alea_raycast_segment_get(result, i, &t_enter, &t_exit, &cell_id, &material_id, &density);
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:d}",
            "t_enter", t_enter, "t_exit", t_exit,
            "cell_id", cell_id, "material_id", material_id,
            "density", density);
        PyList_SET_ITEM(segments, i, item);
    }
    alea_raycast_result_free(result);
    return segments;
}

static PyObject* AleaTHORSystem_ray_first_cell(AleaTHORSystemObject* self, PyObject* args) {
    double ox, oy, oz, dx, dy, dz, t_max = 0.0;
    if (!PyArg_ParseTuple(args, "dddddd|d", &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double t;
    int cell_id = alea_ray_first_cell(self->sys, ox, oy, oz, dx, dy, dz, t_max, &t);
    if (cell_id < 0) Py_RETURN_NONE;
    return Py_BuildValue("(id)", cell_id, t);
}

static PyObject* AleaTHORSystem_find_overlaps(AleaTHORSystemObject* self, PyObject* args) {
    size_t max_pairs = 100;

    if (!PyArg_ParseTuple(args, "|n", &max_pairs)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int* pairs = malloc(max_pairs * 2 * sizeof(int));
    if (!pairs) {
        PyErr_NoMemory();
        return NULL;
    }

    int n;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    n = alea_find_overlaps(self->sys, pairs, max_pairs);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(pairs);
        return NULL;
    }

    PyObject* list = PyList_New(n);
    if (!list) {
        free(pairs);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        PyObject* pair = Py_BuildValue("(ii)", pairs[i*2], pairs[i*2+1]);
        if (!pair) {
            free(pairs);
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, pair);
    }

    free(pairs);
    return list;
}

/* ============================================================================
 * AleaTHORSystem Methods - Volume Estimation
 * ============================================================================ */

static PyObject* AleaTHORSystem_compute_bounding_sphere(AleaTHORSystemObject* self, PyObject* args) {
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "|d", &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double cx, cy, cz, radius;
    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_compute_bounding_sphere(self->sys, tol, &cx, &cy, &cz, &radius);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, "No bounded cells found");
        return NULL;
    }
    return Py_BuildValue("(dddd)", cx, cy, cz, radius);
}

static PyObject* AleaTHORSystem_estimate_cell_volumes(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, radius;
    int n_rays;
    static char* kwlist[] = {"ox", "oy", "oz", "radius", "n_rays", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddi", kwlist,
                                     &ox, &oy, &oz, &radius, &n_rays)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_cell_count(self->sys);
    double* volumes = calloc(count, sizeof(double));
    double* rel_errors = calloc(count, sizeof(double));
    if (!volumes || !rel_errors) {
        free(volumes); free(rel_errors);
        return PyErr_NoMemory();
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_estimate_cell_volumes(self->sys, ox, oy, oz, radius, n_rays,
                                            volumes, rel_errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(volumes); free(rel_errors);
        return NULL;
    }

    if (result < 0) {
        free(volumes); free(rel_errors);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* vol_list = PyList_New(count);
    PyObject* err_list = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(vol_list, i, PyFloat_FromDouble(volumes[i]));
        PyList_SET_ITEM(err_list, i, PyFloat_FromDouble(rel_errors[i]));
    }
    free(volumes);
    free(rel_errors);

    return Py_BuildValue("{s:N, s:N}", "volumes", vol_list, "rel_errors", err_list);
}

static PyObject* AleaTHORSystem_estimate_instance_volumes(AleaTHORSystemObject* self, PyObject* args) {
    int n_rays;
    if (!PyArg_ParseTuple(args, "i", &n_rays)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    size_t count = alea_spatial_index_instance_count(self->sys);
    if (count == 0) {
        PyErr_SetString(PyExc_RuntimeError, "No spatial index built or no instances");
        return NULL;
    }

    double* volumes = calloc(count, sizeof(double));
    double* rel_errors = calloc(count, sizeof(double));
    if (!volumes || !rel_errors) {
        free(volumes); free(rel_errors);
        return PyErr_NoMemory();
    }

    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_estimate_instance_volumes(self->sys, n_rays, volumes, rel_errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(volumes); free(rel_errors);
        return NULL;
    }

    if (result < 0) {
        free(volumes); free(rel_errors);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* vol_list = PyList_New(count);
    PyObject* err_list = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(vol_list, i, PyFloat_FromDouble(volumes[i]));
        PyList_SET_ITEM(err_list, i, PyFloat_FromDouble(rel_errors[i]));
    }
    free(volumes);
    free(rel_errors);

    return Py_BuildValue("{s:N, s:N}", "volumes", vol_list, "rel_errors", err_list);
}

static PyObject* AleaTHORSystem_remove_cells_by_volume(AleaTHORSystemObject* self, PyObject* args) {
    PyObject* vol_list;
    double threshold;
    if (!PyArg_ParseTuple(args, "O!d", &PyList_Type, &vol_list, &threshold)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    Py_ssize_t count = PyList_Size(vol_list);
    double* volumes = malloc(count * sizeof(double));
    if (!volumes) return PyErr_NoMemory();

    for (Py_ssize_t i = 0; i < count; i++) {
        volumes[i] = PyFloat_AsDouble(PyList_GetItem(vol_list, i));
    }
    if (PyErr_Occurred()) { free(volumes); return NULL; }

    int result = alea_remove_cells_by_volume(self->sys, volumes, threshold);
    free(volumes);

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

/* ============================================================================
 * AleaTHORSystem Methods - Cell-Aware Raycast
 * ============================================================================ */

static PyObject* AleaTHORSystem_raycast_cell_aware(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz;
    double t_max = 0.0;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|d", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_raycast_result_t* result = alea_raycast_result_create();
    if (!result) { PyErr_SetString(PyExc_MemoryError, "Failed to allocate raycast result"); return NULL; }

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_raycast_cell_aware(self->sys, ox, oy, oz, dx, dy, dz, t_max, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) { alea_raycast_result_free(result); return NULL; }

    if (rc < 0) {
        alea_raycast_result_free(result);
        PyErr_SetString(PyExc_RuntimeError, "Cell-aware raycast failed");
        return NULL;
    }

    size_t count = alea_raycast_segment_count(result);
    PyObject* segments = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id;
        alea_raycast_segment_get(result, i, &t_enter, &t_exit, &cell_id, &material_id, &density);
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:d}",
            "t_enter", t_enter, "t_exit", t_exit,
            "cell_id", cell_id, "material_id", material_id,
            "density", density);
        PyList_SET_ITEM(segments, i, item);
    }
    alea_raycast_result_free(result);
    return segments;
}
