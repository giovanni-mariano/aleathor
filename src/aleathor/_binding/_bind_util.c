// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Utilities (set_verbose, validate, print_summary, set_tolerance,
 *           clone, reset), config (get_config, set_config, set_log_level).
 */

/* ============================================================================
 * Utilities
 * ============================================================================ */

static PyObject* AleaTHORSystem_set_verbose(AleaTHORSystemObject* self, PyObject* args) {
    int verbose;

    if (!PyArg_ParseTuple(args, "p", &verbose)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.log_level = verbose ? ALEA_LOG_LEVEL_INFO : ALEA_LOG_LEVEL_WARN;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_validate(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int issues = alea_validate(self->sys);
    return PyLong_FromLong(issues);
}

static PyObject* AleaTHORSystem_print_summary(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_print_summary(self->sys);
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_set_tolerance(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double abs_tol = 1e-6;
    double rel_tol = 1e-9;
    double zero_thresh = 1e-10;
    static char* kwlist[] = {"abs_tol", "rel_tol", "zero_thresh", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ddd", kwlist, &abs_tol, &rel_tol, &zero_thresh)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.abs_tol = abs_tol;
    cfg.rel_tol = rel_tol;
    cfg.zero_threshold = zero_thresh;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_clone(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_system_t* cloned = alea_clone(self->sys);
    if (!cloned) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to clone system");
        return NULL;
    }

    AleaTHORSystemObject* obj = (AleaTHORSystemObject*)AleaTHORSystemType.tp_alloc(&AleaTHORSystemType, 0);
    if (!obj) {
        alea_destroy(cloned);
        return NULL;
    }

    obj->sys = cloned;
    obj->owns_sys = 1;
    return (PyObject*)obj;
}

static PyObject* AleaTHORSystem_reset(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_reset(self->sys);
    Py_RETURN_NONE;
}

/* ============================================================================
 * Config
 * ============================================================================ */

static PyObject* AleaTHORSystem_set_log_level(AleaTHORSystemObject* self, PyObject* args) {
    int level;
    if (!PyArg_ParseTuple(args, "i", &level)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);
    cfg.log_level = level;
    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_get_config(AleaTHORSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);

    PyObject* dict = PyDict_New();
    PyDict_SetItemString(dict, "abs_tol", PyFloat_FromDouble(cfg.abs_tol));
    PyDict_SetItemString(dict, "rel_tol", PyFloat_FromDouble(cfg.rel_tol));
    PyDict_SetItemString(dict, "zero_threshold", PyFloat_FromDouble(cfg.zero_threshold));
    PyDict_SetItemString(dict, "dedup", PyBool_FromLong(cfg.dedup));
    PyDict_SetItemString(dict, "log_level", PyLong_FromLong(cfg.log_level));
    PyDict_SetItemString(dict, "surface_policy", PyLong_FromLong(cfg.surface_policy));
    PyDict_SetItemString(dict, "export_materials", PyBool_FromLong(cfg.export_materials));
    PyDict_SetItemString(dict, "export_transforms", PyBool_FromLong(cfg.export_transforms));
    PyDict_SetItemString(dict, "universe_depth", PyLong_FromLong(cfg.universe_depth));
    PyDict_SetItemString(dict, "fill_depth", PyLong_FromLong(cfg.fill_depth));
    PyDict_SetItemString(dict, "trcl_mode", PyLong_FromLong(cfg.trcl_mode));
    PyDict_SetItemString(dict, "transform_mode", PyLong_FromLong(cfg.transform_mode));
    PyDict_SetItemString(dict, "mcnp_max_col", PyLong_FromLong(cfg.mcnp_max_col));
    PyDict_SetItemString(dict, "mcnp_cont_indent", PyLong_FromLong(cfg.mcnp_cont_indent));
    PyDict_SetItemString(dict, "void_max_depth", PyLong_FromLong(cfg.void_max_depth));
    PyDict_SetItemString(dict, "void_min_size", PyFloat_FromDouble(cfg.void_min_size));
    PyDict_SetItemString(dict, "void_samples", PyLong_FromLong(cfg.void_samples));
    PyDict_SetItemString(dict, "merge_cell_weight", PyFloat_FromDouble(cfg.merge_cell_weight));
    PyDict_SetItemString(dict, "merge_surface_weight", PyFloat_FromDouble(cfg.merge_surface_weight));
    PyDict_SetItemString(dict, "merge_max_surfaces", PyLong_FromLong(cfg.merge_max_surfaces));
    PyDict_SetItemString(dict, "merge_min_cells", PyLong_FromLong(cfg.merge_min_cells));
    PyDict_SetItemString(dict, "flatten_max_depth", PyLong_FromLong(cfg.flatten_max_depth));
    return dict;
}

static PyObject* AleaTHORSystem_set_config(AleaTHORSystemObject* self, PyObject* args) {
    PyObject* dict;
    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    alea_config_t cfg = alea_get_config(self->sys);
    PyObject* val;

#define SET_INT(field)    if ((val = PyDict_GetItemString(dict, #field))) cfg.field = (int)PyLong_AsLong(val)
#define SET_DOUBLE(field) if ((val = PyDict_GetItemString(dict, #field))) cfg.field = PyFloat_AsDouble(val)
#define SET_BOOL(field)   if ((val = PyDict_GetItemString(dict, #field))) cfg.field = PyObject_IsTrue(val)

    SET_DOUBLE(abs_tol);
    SET_DOUBLE(rel_tol);
    SET_DOUBLE(zero_threshold);
    SET_BOOL(dedup);
    SET_INT(log_level);
    SET_INT(surface_policy);
    SET_BOOL(export_materials);
    SET_BOOL(export_transforms);
    SET_INT(universe_depth);
    SET_INT(fill_depth);
    SET_INT(trcl_mode);
    SET_INT(transform_mode);
    SET_INT(mcnp_max_col);
    SET_INT(mcnp_cont_indent);
    SET_INT(void_max_depth);
    SET_DOUBLE(void_min_size);
    SET_INT(void_samples);
    SET_DOUBLE(merge_cell_weight);
    SET_DOUBLE(merge_surface_weight);
    SET_INT(merge_max_surfaces);
    SET_INT(merge_min_cells);
    SET_INT(flatten_max_depth);

#undef SET_INT
#undef SET_DOUBLE
#undef SET_BOOL

    if (PyErr_Occurred()) return NULL;

    alea_set_config(self->sys, &cfg);
    Py_RETURN_NONE;
}
