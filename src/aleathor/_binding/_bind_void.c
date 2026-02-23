// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: VoidResult Python type (struct, dealloc, methods, type spec),
 *           generate_void module function.
 */

/* ============================================================================
 * VoidResult Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    void_result_t* result;
    AleaTHORSystemObject* sys_ref;  /* Keep system alive */
} AleaTHORVoidResultObject;

static void AleaTHORVoidResult_dealloc(AleaTHORVoidResultObject* self) {
    if (self->result) {
        alea_void_free(self->result);
    }
    Py_XDECREF(self->sys_ref);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORVoidResult_get_box_count(AleaTHORVoidResultObject* self, void* closure) {
    (void)closure;
    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_void_count(self->result));
}

static PyObject* AleaTHORVoidResult_get_box(AleaTHORVoidResultObject* self, PyObject* args) {
    size_t index;
    if (!PyArg_ParseTuple(args, "n", &index)) return NULL;

    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }

    alea_bbox_t box;
    if (alea_void_get(self->result, index, &box) < 0) {
        PyErr_Format(PyExc_IndexError, "Box index %zu out of range", index);
        return NULL;
    }

    return Py_BuildValue("(dddddd)",
        box.min_x, box.max_x,
        box.min_y, box.max_y,
        box.min_z, box.max_z);
}

static PyObject* AleaTHORVoidResult_get_boxes(AleaTHORVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result) {
        PyErr_SetString(PyExc_RuntimeError, "Result not initialized");
        return NULL;
    }

    size_t count = alea_void_count(self->result);
    PyObject* list = PyList_New(count);
    if (!list) return NULL;

    for (size_t i = 0; i < count; i++) {
        alea_bbox_t box;
        if (alea_void_get(self->result, i, &box) < 0) {
            Py_DECREF(list);
            PyErr_SetString(PyExc_RuntimeError, "Failed to get box");
            return NULL;
        }

        PyObject* bbox = Py_BuildValue("(dddddd)",
            box.min_x, box.max_x,
            box.min_y, box.max_y,
            box.min_z, box.max_z);
        if (!bbox) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, bbox);
    }

    return list;
}

static PyObject* AleaTHORVoidResult_to_node(AleaTHORVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    alea_node_id_t node = alea_void_to_node(self->sys_ref->sys, self->result);
    if (node == ALEA_NODE_ID_INVALID) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create void node");
        return NULL;
    }

    return PyLong_FromUnsignedLong(node);
}

static PyObject* AleaTHORVoidResult_merge(AleaTHORVoidResultObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->result || !self->sys_ref || !self->sys_ref->sys) {
        PyErr_SetString(PyExc_RuntimeError, "Result or system not initialized");
        return NULL;
    }

    int result = alea_void_merge(self->sys_ref->sys, self->result);
    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return PyLong_FromLong(result);
}

static PyGetSetDef AleaTHORVoidResult_getsetters[] = {
    {"box_count", (getter)AleaTHORVoidResult_get_box_count, NULL, "Number of void boxes", NULL},
    {NULL}
};

static PyMethodDef AleaTHORVoidResult_methods[] = {
    {"get_box", (PyCFunction)AleaTHORVoidResult_get_box, METH_VARARGS,
     "get_box(index) -> (xmin, xmax, ymin, ymax, zmin, zmax)"},
    {"get_boxes", (PyCFunction)AleaTHORVoidResult_get_boxes, METH_NOARGS,
     "get_boxes() -> list of bbox tuples"},
    {"to_node", (PyCFunction)AleaTHORVoidResult_to_node, METH_NOARGS,
     "to_node() -> node_id\n\nConvert void result to CSG node (union of boxes)."},
    {"merge", (PyCFunction)AleaTHORVoidResult_merge, METH_NOARGS,
     "merge() -> int\n\nMerge void cells to reduce count while balancing complexity."},
    {NULL}
};

static PyTypeObject AleaTHORVoidResultType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.VoidResult",
    .tp_doc = PyDoc_STR("Result of void generation."),
    .tp_basicsize = sizeof(AleaTHORVoidResultObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)AleaTHORVoidResult_dealloc,
    .tp_methods = AleaTHORVoidResult_methods,
    .tp_getset = AleaTHORVoidResult_getsetters,
};

static PyObject* mod_generate_void(PyObject* self, PyObject* args, PyObject* kwds) {
    (void)self;
    AleaTHORSystemObject* sys_obj;
    PyObject* bounds_obj = Py_None;
    int max_depth = 8;
    double min_size = 0.1;
    int samples_per_node = 27;
    static char* kwlist[] = {"system", "bounds", "max_depth", "min_size", "samples_per_node", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|Oidi", kwlist,
                                     &AleaTHORSystemType, &sys_obj,
                                     &bounds_obj, &max_depth, &min_size, &samples_per_node)) {
        return NULL;
    }

    if (!sys_obj->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Parse bounds if provided */
    alea_bbox_t bounds;
    alea_bbox_t* bounds_ptr = NULL;

    if (bounds_obj != Py_None) {
        if (!PyTuple_Check(bounds_obj) || PyTuple_Size(bounds_obj) != 6) {
            PyErr_SetString(PyExc_TypeError, "bounds must be (xmin, xmax, ymin, ymax, zmin, zmax)");
            return NULL;
        }
        bounds.min_x = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 0));
        bounds.max_x = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 1));
        bounds.min_y = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 2));
        bounds.max_y = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 3));
        bounds.min_z = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 4));
        bounds.max_z = PyFloat_AsDouble(PyTuple_GetItem(bounds_obj, 5));
        if (PyErr_Occurred()) return NULL;
        bounds_ptr = &bounds;
    }

    /* Configure via system config */
    alea_config_t orig_cfg = alea_get_config(sys_obj->sys);
    alea_config_t cfg = orig_cfg;
    cfg.void_max_depth = max_depth;
    cfg.void_min_size = min_size;
    cfg.void_samples = samples_per_node;
    alea_set_config(sys_obj->sys, &cfg);

    void_result_t* result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_void_generate(sys_obj->sys, bounds_ptr);
    Py_END_ALLOW_THREADS

    /* Restore original config */
    alea_set_config(sys_obj->sys, &orig_cfg);

    if (restore_sigint(old_sigint)) {
        if (result) alea_void_free(result);
        return NULL;
    }

    if (!result) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to generate void");
        return NULL;
    }

    /* Create Python object */
    AleaTHORVoidResultObject* obj = PyObject_New(AleaTHORVoidResultObject, &AleaTHORVoidResultType);
    if (!obj) {
        alea_void_free(result);
        return NULL;
    }

    obj->result = result;
    obj->sys_ref = sys_obj;
    Py_INCREF(sys_obj);

    return (PyObject*)obj;
}
