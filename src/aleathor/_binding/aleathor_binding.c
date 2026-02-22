// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file _alea_module.c
 * @brief Python C extension for AleaTHOR-CSG library
 *
 * Uses ONLY the public API from alea.h (alea_* functions).
 *
 * Build: python setup.py build_ext --inplace
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

/* Public API headers */
#include "alea.h"
#include "alea_slice.h"
#include "alea_raycast.h"
#include "alea_mesh.h"

#include <signal.h>

/* alea_log_set_callback is not exposed in the public header (alea.h) but is
 * available in the compiled library.  We forward-declare it here so the
 * Python logging bridge can route C library log messages to Python. */
typedef void (*alea_log_callback_t)(alea_log_level_t level, const char* file,
                                    int line, const char* message, void* user_data);
extern void alea_log_set_callback(alea_log_callback_t callback, void* user_data);

/* SIGINT cooperative-interruption helpers.
 * While the GIL is released Python cannot run its own SIGINT handler.
 * We temporarily install a handler that sets the C library's interrupt flag
 * so that long-running C operations return early with ALEA_ERR_INTERRUPTED. */
static void sigint_handler(int sig) { (void)sig; alea_interrupt(); }

typedef void (*sighandler_func)(int);

static inline sighandler_func install_sigint(void) {
    return signal(SIGINT, sigint_handler);
}

/* Restore the previous handler, check whether the interrupt flag was set,
 * and if so raise KeyboardInterrupt.  Returns 1 if interrupted, 0 otherwise. */
static inline int restore_sigint(sighandler_func old) {
    signal(SIGINT, old);
    int was_interrupted = alea_interrupted();
    alea_clear_interrupt();
    if (was_interrupted) {
        PyErr_SetNone(PyExc_KeyboardInterrupt);
    }
    return was_interrupted;
}

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static PyTypeObject AleaTHORSystemType;
static PyTypeObject AleaTHORVoidResultType;

/* ============================================================================
 * AleaTHORSystem Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_system_t* sys;
    int owns_sys;
} AleaTHORSystemObject;

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
 * AleaTHORSystem Properties
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

static PyGetSetDef AleaTHORSystem_getsetters[] = {
    {"cell_count", (getter)AleaTHORSystem_get_cell_count, NULL, "Number of cells", NULL},
    {"surface_count", (getter)AleaTHORSystem_get_surface_count, NULL, "Number of surfaces", NULL},
    {"universe_count", (getter)AleaTHORSystem_get_universe_count, NULL, "Number of universes", NULL},
    {"spatial_index_instance_count", (getter)AleaTHORSystem_get_spatial_index_instance_count, NULL, "Number of cell instances in spatial index", NULL},
    {NULL}
};

/* ============================================================================
 * AleaTHORSystem Methods - Queries
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
 * AleaTHORSystem Methods - Cell Operations
 * ============================================================================ */

/* Helper: add lattice fields to an existing cell dict.
 * Only adds keys when lat_type != 0 to avoid clutter for non-lattice cells.
 * Returns 0 on success, -1 on failure (with Python exception set). */
static int add_lattice_fields(PyObject* dict, const alea_cell_info_t* info) {
    if (info->lat_type == 0)
        return 0;

    if (PyDict_SetItemString(dict, "lat_type",
            PyLong_FromLong(info->lat_type)) < 0)
        return -1;

    PyObject* dims = Py_BuildValue("(iiiiii)",
        info->lat_fill_dims[0], info->lat_fill_dims[1],
        info->lat_fill_dims[2], info->lat_fill_dims[3],
        info->lat_fill_dims[4], info->lat_fill_dims[5]);
    if (!dims) return -1;
    if (PyDict_SetItemString(dict, "lat_fill_dims", dims) < 0) {
        Py_DECREF(dims); return -1;
    }
    Py_DECREF(dims);

    PyObject* pitch = Py_BuildValue("(ddd)",
        info->lat_pitch[0], info->lat_pitch[1], info->lat_pitch[2]);
    if (!pitch) return -1;
    if (PyDict_SetItemString(dict, "lat_pitch", pitch) < 0) {
        Py_DECREF(pitch); return -1;
    }
    Py_DECREF(pitch);

    PyObject* ll = Py_BuildValue("(ddd)",
        info->lat_lower_left[0], info->lat_lower_left[1], info->lat_lower_left[2]);
    if (!ll) return -1;
    if (PyDict_SetItemString(dict, "lat_lower_left", ll) < 0) {
        Py_DECREF(ll); return -1;
    }
    Py_DECREF(ll);

    PyObject* fill_count = PyLong_FromSize_t(info->lat_fill_count);
    if (!fill_count) return -1;
    if (PyDict_SetItemString(dict, "lat_fill_count", fill_count) < 0) {
        Py_DECREF(fill_count); return -1;
    }
    Py_DECREF(fill_count);

    if (info->lat_fill && info->lat_fill_count > 0) {
        PyObject* fill = PyList_New((Py_ssize_t)info->lat_fill_count);
        if (!fill) return -1;
        for (size_t i = 0; i < info->lat_fill_count; i++) {
            PyList_SET_ITEM(fill, (Py_ssize_t)i, PyLong_FromLong(info->lat_fill[i]));
        }
        if (PyDict_SetItemString(dict, "lat_fill", fill) < 0) {
            Py_DECREF(fill); return -1;
        }
        Py_DECREF(fill);
    } else {
        if (PyDict_SetItemString(dict, "lat_fill", Py_None) < 0)
            return -1;
    }

    return 0;
}

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
 * AleaTHORSystem Methods - Universe Operations
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
 * AleaTHORSystem Methods - Utilities
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
 * Slice Curves API (for matplotlib integration)
 * ============================================================================ */

static const char* curve_type_to_string(alea_curve_type_t type) {
    switch (type) {
        case ALEA_CURVE_LINE: return "line";
        case ALEA_CURVE_LINE_SEGMENT: return "line_segment";
        case ALEA_CURVE_CIRCLE: return "circle";
        case ALEA_CURVE_ARC: return "arc";
        case ALEA_CURVE_ELLIPSE: return "ellipse";
        case ALEA_CURVE_ELLIPSE_ARC: return "ellipse_arc";
        case ALEA_CURVE_POLYGON: return "polygon";
        case ALEA_CURVE_PARALLEL_LINES: return "parallel_lines";
        default: return "none";
    }
}

static PyObject* AleaTHORSystem_get_slice_curves_z(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double z, x_min, x_max, y_min, y_max;

    static char* kwlist[] = {"z", "x_min", "x_max", "y_min", "y_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &z, &x_min, &x_max, &y_min, &y_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, z, x_min, x_max, y_min, y_max);
    alea_slice_curves_t* curves = alea_get_slice_curves(self->sys, &view);
    if (!curves) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) {
        alea_slice_curves_free(curves);
        return NULL;
    }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        PyDict_SetItemString(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        PyDict_SetItemString(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                PyDict_SetItemString(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                if (c.type == ALEA_CURVE_LINE_SEGMENT) {
                    PyDict_SetItemString(dict, "t_min", PyFloat_FromDouble(c.t_min));
                    PyDict_SetItemString(dict, "t_max", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                if (c.type == ALEA_CURVE_ARC) {
                    PyDict_SetItemString(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    PyDict_SetItemString(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                PyDict_SetItemString(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                PyDict_SetItemString(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                PyDict_SetItemString(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                if (c.type == ALEA_CURVE_ELLIPSE_ARC) {
                    PyDict_SetItemString(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    PyDict_SetItemString(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                PyDict_SetItemString(dict, "vertices", verts);
                PyDict_SetItemString(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }

            case ALEA_CURVE_PARALLEL_LINES:
                PyDict_SetItemString(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                PyDict_SetItemString(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;

            default:
                break;
        }

        PyList_SET_ITEM(list, i, dict);
    }

    /* Build result with curves and bounds */
    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "curves", list);
    PyDict_SetItemString(result, "u_min", PyFloat_FromDouble(u_min));
    PyDict_SetItemString(result, "u_max", PyFloat_FromDouble(u_max));
    PyDict_SetItemString(result, "v_min", PyFloat_FromDouble(v_min));
    PyDict_SetItemString(result, "v_max", PyFloat_FromDouble(v_max));
    PyDict_SetItemString(result, "x_min", PyFloat_FromDouble(x_min));
    PyDict_SetItemString(result, "x_max", PyFloat_FromDouble(x_max));
    PyDict_SetItemString(result, "y_min", PyFloat_FromDouble(y_min));
    PyDict_SetItemString(result, "y_max", PyFloat_FromDouble(y_max));

    alea_slice_curves_free(curves);
    return result;
}

static PyObject* AleaTHORSystem_get_slice_curves_y(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double y, x_min, x_max, z_min, z_max;

    static char* kwlist[] = {"y", "x_min", "x_max", "z_min", "z_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &y, &x_min, &x_max, &z_min, &z_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 1, y, x_min, x_max, z_min, z_max);
    alea_slice_curves_t* curves = alea_get_slice_curves(self->sys, &view);
    if (!curves) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Same processing as Z slice */
    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        PyDict_SetItemString(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        PyDict_SetItemString(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                PyDict_SetItemString(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                PyDict_SetItemString(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                PyDict_SetItemString(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                PyDict_SetItemString(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                PyDict_SetItemString(dict, "vertices", verts);
                PyDict_SetItemString(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                PyDict_SetItemString(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                PyDict_SetItemString(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "curves", list);
    PyDict_SetItemString(result, "u_min", PyFloat_FromDouble(u_min));
    PyDict_SetItemString(result, "u_max", PyFloat_FromDouble(u_max));
    PyDict_SetItemString(result, "v_min", PyFloat_FromDouble(v_min));
    PyDict_SetItemString(result, "v_max", PyFloat_FromDouble(v_max));
    /* Viewport bounds for XZ plane (Y slice) */
    PyDict_SetItemString(result, "x_min", PyFloat_FromDouble(x_min));
    PyDict_SetItemString(result, "x_max", PyFloat_FromDouble(x_max));
    PyDict_SetItemString(result, "z_min", PyFloat_FromDouble(z_min));
    PyDict_SetItemString(result, "z_max", PyFloat_FromDouble(z_max));

    alea_slice_curves_free(curves);
    return result;
}

static PyObject* AleaTHORSystem_get_slice_curves_x(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double x, y_min, y_max, z_min, z_max;

    static char* kwlist[] = {"x", "y_min", "y_max", "z_min", "z_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &x, &y_min, &y_max, &z_min, &z_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 0, x, y_min, y_max, z_min, z_max);
    alea_slice_curves_t* curves = alea_get_slice_curves(self->sys, &view);
    if (!curves) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        PyDict_SetItemString(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        PyDict_SetItemString(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                PyDict_SetItemString(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                PyDict_SetItemString(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                PyDict_SetItemString(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                PyDict_SetItemString(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                PyDict_SetItemString(dict, "vertices", verts);
                PyDict_SetItemString(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                PyDict_SetItemString(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                PyDict_SetItemString(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "curves", list);
    PyDict_SetItemString(result, "u_min", PyFloat_FromDouble(u_min));
    PyDict_SetItemString(result, "u_max", PyFloat_FromDouble(u_max));
    PyDict_SetItemString(result, "v_min", PyFloat_FromDouble(v_min));
    PyDict_SetItemString(result, "v_max", PyFloat_FromDouble(v_max));
    /* Viewport bounds for YZ plane (X slice) */
    PyDict_SetItemString(result, "y_min", PyFloat_FromDouble(y_min));
    PyDict_SetItemString(result, "y_max", PyFloat_FromDouble(y_max));
    PyDict_SetItemString(result, "z_min", PyFloat_FromDouble(z_min));
    PyDict_SetItemString(result, "z_max", PyFloat_FromDouble(z_max));

    alea_slice_curves_free(curves);
    return result;
}

/* Arbitrary plane slice curves */

static PyObject* AleaTHORSystem_get_slice_curves(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj;
    PyObject* normal_obj;
    PyObject* up_obj;
    double u_min, u_max, v_min, v_max;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdddd", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Parse tuples */
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) {
        PyErr_SetString(PyExc_TypeError, "origin must be (x, y, z) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) {
        PyErr_SetString(PyExc_TypeError, "normal must be (nx, ny, nz) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError, "up must be (ux, uy, uz) tuple");
        return NULL;
    }

    /* Initialize view */
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    /* Get curves */
    alea_slice_curves_t* curves = alea_get_slice_curves(self->sys, &view);
    if (!curves) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Build result (same as axis-aligned versions) */
    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        PyDict_SetItemString(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        PyDict_SetItemString(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                PyDict_SetItemString(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                if (c.type == ALEA_CURVE_ARC) {
                    PyDict_SetItemString(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    PyDict_SetItemString(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                PyDict_SetItemString(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                PyDict_SetItemString(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                PyDict_SetItemString(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                PyDict_SetItemString(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                if (c.type == ALEA_CURVE_ELLIPSE_ARC) {
                    PyDict_SetItemString(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    PyDict_SetItemString(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                PyDict_SetItemString(dict, "vertices", verts);
                PyDict_SetItemString(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                PyDict_SetItemString(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                PyDict_SetItemString(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                PyDict_SetItemString(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double cu_min, cu_max, cv_min, cv_max;
    alea_slice_curves_bounds(curves, &cu_min, &cu_max, &cv_min, &cv_max);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "curves", list);
    PyDict_SetItemString(result, "u_min", PyFloat_FromDouble(cu_min));
    PyDict_SetItemString(result, "u_max", PyFloat_FromDouble(cu_max));
    PyDict_SetItemString(result, "v_min", PyFloat_FromDouble(cv_min));
    PyDict_SetItemString(result, "v_max", PyFloat_FromDouble(cv_max));
    /* Viewport bounds */
    PyDict_SetItemString(result, "view_u_min", PyFloat_FromDouble(u_min));
    PyDict_SetItemString(result, "view_u_max", PyFloat_FromDouble(u_max));
    PyDict_SetItemString(result, "view_v_min", PyFloat_FromDouble(v_min));
    PyDict_SetItemString(result, "view_v_max", PyFloat_FromDouble(v_max));

    alea_slice_curves_free(curves);
    return result;
}

/* Arbitrary plane grid query */

static PyObject* AleaTHORSystem_find_cells_grid(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj;
    PyObject* normal_obj;
    PyObject* up_obj;
    double u_min, u_max, v_min, v_max;
    int nu, nv;
    int universe_depth = -1;  /* Default: innermost cell */
    int detect_errors = 0;    /* Default: no error detection */

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "nu", "nv",
                             "universe_depth", "detect_errors", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddii|ip", kwlist,
            &origin_obj, &normal_obj, &up_obj,
            &u_min, &u_max, &v_min, &v_max, &nu, &nv,
            &universe_depth, &detect_errors)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Parse tuples */
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) {
        PyErr_SetString(PyExc_TypeError, "origin must be (x, y, z) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) {
        PyErr_SetString(PyExc_TypeError, "normal must be (nx, ny, nz) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError, "up must be (ux, uy, uz) tuple");
        return NULL;
    }

    if (nu <= 0 || nv <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    /* Initialize view */
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    /* Allocate arrays */
    int* cell_ids = malloc(nu * nv * sizeof(int));
    int* mat_ids = malloc(nu * nv * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nu * nv * sizeof(uint8_t)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && !errors)) {
        free(cell_ids); free(mat_ids); free(errors);
        return PyErr_NoMemory();
    }

    if (alea_find_cells_grid(self->sys, &view, nu, nv,
                                     universe_depth, cell_ids, mat_ids, errors) < 0) {
        free(cell_ids); free(mat_ids); free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    /* Build result lists */
    PyObject* cell_list = PyList_New(nu * nv);
    PyObject* mat_list = PyList_New(nu * nv);
    for (int i = 0; i < nu * nv; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = PyList_New(nu * nv);
        for (int i = 0; i < nu * nv; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }

    free(cell_ids);
    free(mat_ids);
    free(errors);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    PyDict_SetItemString(result, "nu", PyLong_FromLong(nu));
    PyDict_SetItemString(result, "nv", PyLong_FromLong(nv));
    PyDict_SetItemString(result, "u_min", PyFloat_FromDouble(u_min));
    PyDict_SetItemString(result, "u_max", PyFloat_FromDouble(u_max));
    PyDict_SetItemString(result, "v_min", PyFloat_FromDouble(v_min));
    PyDict_SetItemString(result, "v_max", PyFloat_FromDouble(v_max));

    return result;
}

/* Grid cell queries for fills */

static PyObject* AleaTHORSystem_find_cells_grid_z(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double z, x_min, x_max, y_min, y_max;
    int nx, ny;
    int universe_depth = -1;  /* Default: innermost cell */
    int detect_errors = 0;    /* Default: no error detection */

    static char* kwlist[] = {"z", "x_min", "x_max", "y_min", "y_max", "nx", "ny",
                             "universe_depth", "detect_errors", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|ip", kwlist,
            &z, &x_min, &x_max, &y_min, &y_max, &nx, &ny,
            &universe_depth, &detect_errors)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (nx <= 0 || ny <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(nx * ny * sizeof(int));
    int* mat_ids = malloc(nx * ny * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nx * ny * sizeof(uint8_t)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && !errors)) {
        free(cell_ids); free(mat_ids); free(errors);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, z, x_min, x_max, y_min, y_max);
    if (alea_find_cells_grid(self->sys, &view, nx, ny,
                                     universe_depth, cell_ids, mat_ids, errors) < 0) {
        free(cell_ids); free(mat_ids); free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    /* Create lists for results */
    PyObject* cell_list = PyList_New(nx * ny);
    PyObject* mat_list = PyList_New(nx * ny);
    for (int i = 0; i < nx * ny; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = PyList_New(nx * ny);
        for (int i = 0; i < nx * ny; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }

    free(cell_ids);
    free(mat_ids);
    free(errors);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    PyDict_SetItemString(result, "nx", PyLong_FromLong(nx));
    PyDict_SetItemString(result, "ny", PyLong_FromLong(ny));
    PyDict_SetItemString(result, "x_min", PyFloat_FromDouble(x_min));
    PyDict_SetItemString(result, "x_max", PyFloat_FromDouble(x_max));
    PyDict_SetItemString(result, "y_min", PyFloat_FromDouble(y_min));
    PyDict_SetItemString(result, "y_max", PyFloat_FromDouble(y_max));

    return result;
}

static PyObject* AleaTHORSystem_find_cells_grid_y(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double y, x_min, x_max, z_min, z_max;
    int nx, nz;
    int universe_depth = -1;  /* Default: innermost cell */
    int detect_errors = 0;    /* Default: no error detection */

    static char* kwlist[] = {"y", "x_min", "x_max", "z_min", "z_max", "nx", "nz",
                             "universe_depth", "detect_errors", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|ip", kwlist,
            &y, &x_min, &x_max, &z_min, &z_max, &nx, &nz,
            &universe_depth, &detect_errors)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (nx <= 0 || nz <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(nx * nz * sizeof(int));
    int* mat_ids = malloc(nx * nz * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nx * nz * sizeof(uint8_t)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && !errors)) {
        free(cell_ids); free(mat_ids); free(errors);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 1, y, x_min, x_max, z_min, z_max);
    if (alea_find_cells_grid(self->sys, &view, nx, nz,
                                     universe_depth, cell_ids, mat_ids, errors) < 0) {
        free(cell_ids); free(mat_ids); free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* cell_list = PyList_New(nx * nz);
    PyObject* mat_list = PyList_New(nx * nz);
    for (int i = 0; i < nx * nz; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = PyList_New(nx * nz);
        for (int i = 0; i < nx * nz; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }

    free(cell_ids);
    free(mat_ids);
    free(errors);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    PyDict_SetItemString(result, "nx", PyLong_FromLong(nx));
    PyDict_SetItemString(result, "nz", PyLong_FromLong(nz));
    PyDict_SetItemString(result, "x_min", PyFloat_FromDouble(x_min));
    PyDict_SetItemString(result, "x_max", PyFloat_FromDouble(x_max));
    PyDict_SetItemString(result, "z_min", PyFloat_FromDouble(z_min));
    PyDict_SetItemString(result, "z_max", PyFloat_FromDouble(z_max));

    return result;
}

static PyObject* AleaTHORSystem_find_cells_grid_x(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    double x, y_min, y_max, z_min, z_max;
    int ny, nz;
    int universe_depth = -1;  /* Default: innermost cell */
    int detect_errors = 0;    /* Default: no error detection */

    static char* kwlist[] = {"x", "y_min", "y_max", "z_min", "z_max", "ny", "nz",
                             "universe_depth", "detect_errors", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|ip", kwlist,
            &x, &y_min, &y_max, &z_min, &z_max, &ny, &nz,
            &universe_depth, &detect_errors)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (ny <= 0 || nz <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(ny * nz * sizeof(int));
    int* mat_ids = malloc(ny * nz * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(ny * nz * sizeof(uint8_t)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && !errors)) {
        free(cell_ids); free(mat_ids); free(errors);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 0, x, y_min, y_max, z_min, z_max);
    if (alea_find_cells_grid(self->sys, &view, ny, nz,
                                     universe_depth, cell_ids, mat_ids, errors) < 0) {
        free(cell_ids); free(mat_ids); free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* cell_list = PyList_New(ny * nz);
    PyObject* mat_list = PyList_New(ny * nz);
    for (int i = 0; i < ny * nz; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = PyList_New(ny * nz);
        for (int i = 0; i < ny * nz; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }

    free(cell_ids);
    free(mat_ids);
    free(errors);

    PyObject* result = PyDict_New();
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    PyDict_SetItemString(result, "ny", PyLong_FromLong(ny));
    PyDict_SetItemString(result, "nz", PyLong_FromLong(nz));
    PyDict_SetItemString(result, "y_min", PyFloat_FromDouble(y_min));
    PyDict_SetItemString(result, "y_max", PyFloat_FromDouble(y_max));
    PyDict_SetItemString(result, "z_min", PyFloat_FromDouble(z_min));
    PyDict_SetItemString(result, "z_max", PyFloat_FromDouble(z_max));

    return result;
}

/* Label positioning functions */

static PyObject* AleaTHORSystem_find_label_positions(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* ids_obj;
    int width, height;
    int min_pixels = 100;

    static char* kwlist[] = {"ids", "width", "height", "min_pixels", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oii|i", kwlist,
            &ids_obj, &width, &height, &min_pixels)) return NULL;

    if (width <= 0 || height <= 0) {
        PyErr_SetString(PyExc_ValueError, "Width and height must be positive");
        return NULL;
    }

    /* Convert Python list to C array */
    if (!PyList_Check(ids_obj)) {
        PyErr_SetString(PyExc_TypeError, "ids must be a list");
        return NULL;
    }

    Py_ssize_t list_size = PyList_Size(ids_obj);
    if (list_size != width * height) {
        PyErr_SetString(PyExc_ValueError, "ids list size must equal width * height");
        return NULL;
    }

    int* ids = malloc(width * height * sizeof(int));
    if (!ids) return PyErr_NoMemory();

    for (Py_ssize_t i = 0; i < list_size; i++) {
        PyObject* item = PyList_GET_ITEM(ids_obj, i);
        ids[i] = (int)PyLong_AsLong(item);
    }

    /* Call C function */
    alea_label_position_t* labels = NULL;
    int count = 0;

    int rc = alea_find_label_positions(ids, width, height, min_pixels, &labels, &count);
    free(ids);

    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Label position computation failed");
        return NULL;
    }

    /* Build result list */
    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyObject* label_dict = PyDict_New();
        PyDict_SetItemString(label_dict, "id", PyLong_FromLong(labels[i].id));
        PyDict_SetItemString(label_dict, "px", PyLong_FromLong(labels[i].px));
        PyDict_SetItemString(label_dict, "py", PyLong_FromLong(labels[i].py));
        PyDict_SetItemString(label_dict, "pixel_count", PyLong_FromLong(labels[i].pixel_count));
        PyList_SET_ITEM(result, i, label_dict);
    }

    free(labels);
    return result;
}


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

/* ============================================================================
 * AleaTHORSystem Methods - Grid Overlap Check
 * ============================================================================ */

static PyObject* AleaTHORSystem_check_grid_overlaps(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    int nu, nv;
    int universe_depth = -1;
    PyObject* cell_ids_obj;
    PyObject* errors_obj;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
                             "nu", "nv", "cell_ids", "errors", "universe_depth", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddiiOO|i", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &nu, &nv, &cell_ids_obj, &errors_obj, &universe_depth)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) return NULL;
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) return NULL;
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    if (!PyList_Check(cell_ids_obj) || !PyList_Check(errors_obj)) {
        PyErr_SetString(PyExc_TypeError, "cell_ids and errors must be lists");
        return NULL;
    }

    Py_ssize_t size = nu * nv;
    int* cell_ids = malloc(size * sizeof(int));
    uint8_t* errors = malloc(size * sizeof(uint8_t));
    if (!cell_ids || !errors) { free(cell_ids); free(errors); return PyErr_NoMemory(); }

    for (Py_ssize_t i = 0; i < size; i++) {
        cell_ids[i] = (int)PyLong_AsLong(PyList_GET_ITEM(cell_ids_obj, i));
        errors[i] = (uint8_t)PyLong_AsLong(PyList_GET_ITEM(errors_obj, i));
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    int rc = alea_check_grid_overlaps(self->sys, &view, nu, nv,
                                           universe_depth, cell_ids, errors);
    free(cell_ids);

    if (rc < 0) {
        free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Overlap check failed");
        return NULL;
    }

    /* Return updated errors list */
    PyObject* result = PyList_New(size);
    for (Py_ssize_t i = 0; i < size; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(errors[i]));
    }
    free(errors);
    return result;
}

/* ============================================================================
 * AleaTHORSystem Methods - Surface Label Positions
 * ============================================================================ */

static PyObject* AleaTHORSystem_find_surface_label_positions(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    int width, height, margin = 20;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
                             "width", "height", "margin", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddii|i", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &width, &height, &margin)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) return NULL;
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) return NULL;
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    alea_slice_curves_t* curves = alea_get_slice_curves(self->sys, &view);
    if (!curves) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc = alea_find_surface_label_positions(curves, u_min, u_max, v_min, v_max,
                                                    width, height, margin, &labels, &count);
    alea_slice_curves_free(curves);

    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, "Surface label position computation failed");
        return NULL;
    }

    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyObject* d = PyDict_New();
        PyDict_SetItemString(d, "id", PyLong_FromLong(labels[i].id));
        PyDict_SetItemString(d, "px", PyLong_FromLong(labels[i].px));
        PyDict_SetItemString(d, "py", PyLong_FromLong(labels[i].py));
        PyList_SET_ITEM(result, i, d);
    }
    free(labels);
    return result;
}

/* ============================================================================
 * AleaTHORSystem Methods - Logging / Debug
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

/* ============================================================================
 * AleaTHORSystem Methods - Mesh Module
 * ============================================================================ */

static PyObject* AleaTHORSystem_mesh_export(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    int nx = 10, ny = 10, nz = 10;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    const char* format_str = "gmsh";
    int void_material_id = 0;
    double auto_pad = 0.01;

    static char* kwlist[] = {"filename", "nx", "ny", "nz",
                             "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
                             "format", "void_material_id", "auto_pad", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|iiiddddddsid", kwlist,
            &filename, &nx, &ny, &nz,
            &x_min, &x_max, &y_min, &y_max, &z_min, &z_max,
            &format_str, &void_material_id, &auto_pad)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = nx; cfg.ny = ny; cfg.nz = nz;
    cfg.x_min = x_min; cfg.x_max = x_max;
    cfg.y_min = y_min; cfg.y_max = y_max;
    cfg.z_min = z_min; cfg.z_max = z_max;
    cfg.void_material_id = void_material_id;
    cfg.auto_pad = auto_pad;

    if (strcmp(format_str, "vtk") == 0) {
        cfg.format = ALEA_MESH_VTK;
    } else {
        cfg.format = ALEA_MESH_GMSH;
    }

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_mesh_export_system(self->sys, &cfg, filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* AleaTHORSystem_mesh_sample(AleaTHORSystemObject* self, PyObject* args, PyObject* kwds) {
    int nx = 10, ny = 10, nz = 10;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    int void_material_id = 0;
    double auto_pad = 0.01;

    static char* kwlist[] = {"nx", "ny", "nz",
                             "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
                             "void_material_id", "auto_pad", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiiddddddid", kwlist,
            &nx, &ny, &nz,
            &x_min, &x_max, &y_min, &y_max, &z_min, &z_max,
            &void_material_id, &auto_pad)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = nx; cfg.ny = ny; cfg.nz = nz;
    cfg.x_min = x_min; cfg.x_max = x_max;
    cfg.y_min = y_min; cfg.y_max = y_max;
    cfg.z_min = z_min; cfg.z_max = z_max;
    cfg.void_material_id = void_material_id;
    cfg.auto_pad = auto_pad;

    alea_mesh_result_t* mesh;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    mesh = alea_mesh_sample(self->sys, &cfg);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (!mesh) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Build result dict */
    Py_ssize_t total = (Py_ssize_t)mesh->nx * mesh->ny * mesh->nz;

    PyObject* mat_ids = PyList_New(total);
    PyObject* cell_ids = PyList_New(total);
    if (!mat_ids || !cell_ids) {
        Py_XDECREF(mat_ids); Py_XDECREF(cell_ids);
        alea_mesh_result_free(mesh);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < total; i++) {
        PyList_SET_ITEM(mat_ids, i, PyLong_FromLong(mesh->material_ids[i]));
        PyList_SET_ITEM(cell_ids, i, PyLong_FromLong(mesh->cell_ids[i]));
    }

    /* Node positions */
    PyObject* x_nodes = PyList_New(mesh->nx + 1);
    PyObject* y_nodes = PyList_New(mesh->ny + 1);
    PyObject* z_nodes = PyList_New(mesh->nz + 1);
    if (!x_nodes || !y_nodes || !z_nodes) {
        Py_XDECREF(mat_ids); Py_XDECREF(cell_ids);
        Py_XDECREF(x_nodes); Py_XDECREF(y_nodes); Py_XDECREF(z_nodes);
        alea_mesh_result_free(mesh);
        return NULL;
    }
    for (int i = 0; i <= mesh->nx; i++)
        PyList_SET_ITEM(x_nodes, i, PyFloat_FromDouble(mesh->x_nodes[i]));
    for (int i = 0; i <= mesh->ny; i++)
        PyList_SET_ITEM(y_nodes, i, PyFloat_FromDouble(mesh->y_nodes[i]));
    for (int i = 0; i <= mesh->nz; i++)
        PyList_SET_ITEM(z_nodes, i, PyFloat_FromDouble(mesh->z_nodes[i]));

    PyObject* result = Py_BuildValue("{s:N, s:N, s:N, s:N, s:N, s:i, s:i, s:i}",
        "material_ids", mat_ids,
        "cell_ids", cell_ids,
        "x_nodes", x_nodes,
        "y_nodes", y_nodes,
        "z_nodes", z_nodes,
        "nx", mesh->nx,
        "ny", mesh->ny,
        "nz", mesh->nz);

    alea_mesh_result_free(mesh);
    return result;
}

/* ============================================================================
 * AleaTHORSystem Method Table
 * ============================================================================ */

// static PyMethodDef AleaTHORSystem_methods[] = {
//     /* Queries */
//     {"find_cell", (PyCFunction)AleaTHORSystem_find_cell, METH_VARARGS,
//      "find_cell(x, y, z) -> (cell_id, material_id) or None\n\nFind cell containing point."},
//     {"point_inside", (PyCFunction)AleaTHORSystem_point_inside, METH_VARARGS,
//      "point_inside(node_id, x, y, z) -> bool\n\nTest if point is inside CSG tree."},
//     {"material_at", (PyCFunction)AleaTHORSystem_material_at, METH_VARARGS,
//      "material_at(x, y, z) -> int or None\n\nGet material ID at point."},
//     {"find_overlaps", (PyCFunction)AleaTHORSystem_find_overlaps, METH_VARARGS,
//      "find_overlaps(max_pairs=100) -> list of (cell_idx, cell_idx)\n\nFind overlapping cells."},

//     /* Cell operations */
//     {"get_cell", (PyCFunction)AleaTHORSystem_get_cell, METH_VARARGS,
//      "get_cell(cell_id) -> dict\n\nGet cell by MCNP ID."},
//     {"get_cell_by_index", (PyCFunction)AleaTHORSystem_get_cell_by_index, METH_VARARGS,
//      "get_cell_by_index(index) -> dict\n\nGet cell by array index."},
//     {"get_cells", (PyCFunction)AleaTHORSystem_get_cells, METH_NOARGS,
//      "get_cells() -> list of dict\n\nGet all cells."},

//     /* Universe operations */
//     {"build_universe_index", (PyCFunction)AleaTHORSystem_build_universe_index, METH_NOARGS,
//      "build_universe_index()\n\nBuild universe lookup tables."},
//     {"flatten_universe", (PyCFunction)AleaTHORSystem_flatten_universe, METH_VARARGS,
//      "flatten_universe(universe_id=0) -> int\n\nFlatten universe hierarchy."},
//     {"get_universe", (PyCFunction)AleaTHORSystem_get_universe, METH_VARARGS,
//      "get_universe(universe_id) -> dict\n\nGet universe info."},

//     /* Primitive creation - sense: -1=inside/negative, +1=outside/positive */
//     {"create_plane", (PyCFunction)AleaTHORSystem_create_plane, METH_VARARGS,
//      "create_plane(a, b, c, d, sense) -> node_id\n\nCreate plane halfspace. sense: -1 or +1"},
//     {"create_sphere", (PyCFunction)AleaTHORSystem_create_sphere, METH_VARARGS,
//      "create_sphere(cx, cy, cz, radius, sense) -> node_id\n\nCreate sphere halfspace. sense: -1=inside, +1=outside"},
//     {"create_box", (PyCFunction)AleaTHORSystem_create_box, METH_VARARGS,
//      "create_box(xmin, xmax, ymin, ymax, zmin, zmax, sense) -> node_id\n\nCreate box halfspace. sense: -1=inside, +1=outside"},
//     {"create_cylinder_z", (PyCFunction)AleaTHORSystem_create_cylinder_z, METH_VARARGS,
//      "create_cylinder_z(cx, cy, radius, sense) -> node_id\n\nCreate Z-cylinder halfspace. sense: -1=inside, +1=outside"},

//     /* Boolean operations */
//     {"create_union", (PyCFunction)AleaTHORSystem_create_union, METH_VARARGS,
//      "create_union(a, b) -> node_id\n\nCreate union of two nodes."},
//     {"create_intersection", (PyCFunction)AleaTHORSystem_create_intersection, METH_VARARGS,
//      "create_intersection(a, b) -> node_id\n\nCreate intersection of two nodes."},
//     {"create_difference", (PyCFunction)AleaTHORSystem_create_difference, METH_VARARGS,
//      "create_difference(a, b) -> node_id\n\nCreate difference (a - b)."},
//     {"create_complement", (PyCFunction)AleaTHORSystem_create_complement, METH_VARARGS,
//      "create_complement(a) -> node_id\n\nCreate complement (not a)."},
//     {"create_union_many", (PyCFunction)AleaTHORSystem_create_union_many, METH_VARARGS,
//      "create_union_many([nodes]) -> node_id\n\nCreate union of multiple nodes."},
//     {"create_intersection_many", (PyCFunction)AleaTHORSystem_create_intersection_many, METH_VARARGS,
//      "create_intersection_many([nodes]) -> node_id\n\nCreate intersection of multiple nodes."},

//     /* Cell registration */
//     {"set_sense", (PyCFunction)AleaTHORSystem_set_sense, METH_VARARGS,
//      "set_sense(node_id, sense)\n\nSet primitive sense (+1 outside, -1 inside)."},
//     {"register_surface", (PyCFunction)AleaTHORSystem_register_surface, METH_VARARGS,
//      "register_surface(node_id, surface_id)\n\nRegister primitive as MCNP surface."},
//     {"register_cell", (PyCFunctionWithKeywords)AleaTHORSystem_register_cell, METH_VARARGS | METH_KEYWORDS,
//      "register_cell(cell_id, root_node, material_id=0, density=0.0, universe_id=0) -> index"},

//     /* Export */
//     {"export_mcnp", (PyCFunction)AleaTHORSystem_export_mcnp, METH_VARARGS | METH_KEYWORDS,
//      "export_mcnp(filename, deduplicate=True)\n\nExport to MCNP format."},
//     {"export_openmc", (PyCFunction)AleaTHORSystem_export_openmc, METH_VARARGS | METH_KEYWORDS,
//      "export_openmc(filename)\n\nExport to OpenMC XML format."},

//     /* Merge */
//     {"merge", (PyCFunction)AleaTHORSystem_merge, METH_VARARGS | METH_KEYWORDS,
//      "merge(other, id_offset=0)\n\nMerge another system into this one."},

//     /* Utilities */
//     {"set_verbose", (PyCFunction)AleaTHORSystem_set_verbose, METH_VARARGS,
//      "set_verbose(enabled)\n\nEnable/disable verbose output."},
//     {"validate", (PyCFunction)AleaTHORSystem_validate, METH_NOARGS,
//      "validate() -> int\n\nValidate system integrity. Returns number of issues."},
//     {"print_summary", (PyCFunction)AleaTHORSystem_print_summary, METH_NOARGS,
//      "print_summary()\n\nPrint system summary to stdout."},
//     {"set_tolerance", (PyCFunction)AleaTHORSystem_set_tolerance, METH_VARARGS | METH_KEYWORDS,
//      "set_tolerance(abs_tol=1e-6, rel_tol=1e-9, zero_thresh=1e-10)\n\nSet deduplication tolerances."},
//     {"clone", (PyCFunction)AleaTHORSystem_clone, METH_NOARGS,
//      "clone() -> System\n\nCreate deep copy of system."},
//     {"reset", (PyCFunction)AleaTHORSystem_reset, METH_NOARGS,
//      "reset()\n\nReset system to empty state."},

//     {NULL}
// };


/* ============================================================================
 * CSG Node Tree Inspection
 * ============================================================================ */

static PyObject* build_node_tree(const alea_system_t* sys, alea_node_id_t node) {
    alea_operation_t op = alea_node_operation(sys, node);

    if (op == ALEA_OP_PRIMITIVE) {
        int sid = alea_node_surface_id(sys, node);
        int sense = alea_node_sense(sys, node);
        return Py_BuildValue("(iii)", (int)op, sid, sense);
    }

    if (op == ALEA_OP_COMPLEMENT) {
        alea_node_id_t child = alea_node_left(sys, node);
        PyObject* child_tree = build_node_tree(sys, child);
        if (!child_tree) return NULL;
        PyObject* result = PyTuple_New(2);
        if (!result) { Py_DECREF(child_tree); return NULL; }
        PyTuple_SET_ITEM(result, 0, PyLong_FromLong((int)op));
        PyTuple_SET_ITEM(result, 1, child_tree);
        return result;
    }

    /* Binary: UNION, INTERSECTION, DIFFERENCE */
    alea_node_id_t left = alea_node_left(sys, node);
    alea_node_id_t right = alea_node_right(sys, node);
    PyObject* left_tree = build_node_tree(sys, left);
    if (!left_tree) return NULL;
    PyObject* right_tree = build_node_tree(sys, right);
    if (!right_tree) { Py_DECREF(left_tree); return NULL; }
    PyObject* result = PyTuple_New(3);
    if (!result) { Py_DECREF(left_tree); Py_DECREF(right_tree); return NULL; }
    PyTuple_SET_ITEM(result, 0, PyLong_FromLong((int)op));
    PyTuple_SET_ITEM(result, 1, left_tree);
    PyTuple_SET_ITEM(result, 2, right_tree);
    return result;
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


static PyMethodDef AleaTHORSystem_methods[] = {
    /* Queries */
    {"find_cell", (PyCFunction)AleaTHORSystem_find_cell, METH_VARARGS,
     "find_cell(x, y, z) -> (cell_id, material_id) or None\n\nFind cell containing point."},
    {"point_inside", (PyCFunction)AleaTHORSystem_point_inside, METH_VARARGS,
     "point_inside(node_id, x, y, z) -> bool\n\nTest if point is inside CSG tree."},
    {"material_at", (PyCFunction)AleaTHORSystem_material_at, METH_VARARGS,
     "material_at(x, y, z) -> int or None\n\nGet material ID at point."},
    {"find_overlaps", (PyCFunction)AleaTHORSystem_find_overlaps, METH_VARARGS,
     "find_overlaps(max_pairs=100) -> list of (cell_idx, cell_idx)\n\nFind overlapping cells."},

    /* Cell operations */
    {"get_cell", (PyCFunction)AleaTHORSystem_get_cell, METH_VARARGS,
     "get_cell(cell_id) -> dict\n\nGet cell by MCNP ID."},
    {"get_cell_by_index", (PyCFunction)AleaTHORSystem_get_cell_by_index, METH_VARARGS,
     "get_cell_by_index(index) -> dict\n\nGet cell by array index."},
    {"get_cells", (PyCFunction)AleaTHORSystem_get_cells, METH_NOARGS,
     "get_cells() -> list of dict\n\nGet all cells."},
    {"cell_find", (PyCFunction)AleaTHORSystem_cell_find, METH_VARARGS,
     "cell_find(cell_id) -> int or None\n\nFind cell index by MCNP cell ID. O(1) lookup."},
    {"find_all_cells", (PyCFunction)AleaTHORSystem_find_all_cells, METH_VARARGS,
     "find_all_cells(x, y, z) -> list of dict\n\n"
     "Find all cells at a point across hierarchy depths.\n"
     "Returns list of dicts with cell_id, cell_index, material_id, universe_id, fill_universe, depth, local_x/y/z."},

    /* CSG Node Inspection */
    {"node_tree", (PyCFunction)AleaTHORSystem_node_tree, METH_VARARGS,
     "node_tree(node_id) -> tuple\n\n"
     "Walk CSG tree and return nested tuple structure.\n"
     "Primitive: (0, surface_id, sense)\n"
     "Union: (1, left, right)\n"
     "Intersection: (2, left, right)\n"
     "Difference: (3, left, right)\n"
     "Complement: (4, child)\n"},

    /* Universe operations */
    {"build_universe_index", (PyCFunction)AleaTHORSystem_build_universe_index, METH_NOARGS,
     "build_universe_index()\n\nBuild universe lookup tables."},
    {"flatten_universe", (PyCFunction)AleaTHORSystem_flatten_universe, METH_VARARGS,
     "flatten_universe(universe_id=0) -> int\n\nFlatten universe hierarchy."},
    {"build_spatial_index", (PyCFunction)AleaTHORSystem_build_spatial_index, METH_NOARGS,
     "build_spatial_index()\n\nBuild spatial index for fast queries (no flattening needed)."},
    {"get_universe", (PyCFunction)AleaTHORSystem_get_universe, METH_VARARGS,
     "get_universe(universe_id) -> dict\n\nGet universe info."},

    /* Primitive creation */
    {"create_plane", (PyCFunction)AleaTHORSystem_create_plane, METH_VARARGS,
     "create_plane(a, b, c, d, sense) -> node_id\n\nCreate plane halfspace."},
    {"create_sphere", (PyCFunction)AleaTHORSystem_create_sphere, METH_VARARGS,
     "create_sphere(cx, cy, cz, radius, sense) -> node_id\n\nCreate sphere halfspace."},
    {"create_box", (PyCFunction)AleaTHORSystem_create_box, METH_VARARGS,
     "create_box(xmin, xmax, ymin, ymax, zmin, zmax, sense) -> node_id\n\nCreate box halfspace."},
    {"create_cylinder_z", (PyCFunction)AleaTHORSystem_create_cylinder_z, METH_VARARGS,
     "create_cylinder_z(cx, cy, radius, sense) -> node_id\n\nCreate Z-cylinder halfspace."},

    /* Boolean operations */
    {"create_union", (PyCFunction)AleaTHORSystem_create_union, METH_VARARGS,
     "create_union(a, b) -> node_id\n\nCreate union of two nodes."},
    {"create_intersection", (PyCFunction)AleaTHORSystem_create_intersection, METH_VARARGS,
     "create_intersection(a, b) -> node_id\n\nCreate intersection of two nodes."},
    {"create_difference", (PyCFunction)AleaTHORSystem_create_difference, METH_VARARGS,
     "create_difference(a, b) -> node_id\n\nCreate difference (a - b)."},
    {"create_complement", (PyCFunction)AleaTHORSystem_create_complement, METH_VARARGS,
     "create_complement(a) -> node_id\n\nCreate complement (not a)."},
    {"create_union_many", (PyCFunction)AleaTHORSystem_create_union_many, METH_VARARGS,
     "create_union_many([nodes]) -> node_id\n\nCreate union of multiple nodes."},
    {"create_intersection_many", (PyCFunction)AleaTHORSystem_create_intersection_many, METH_VARARGS,
     "create_intersection_many([nodes]) -> node_id\n\nCreate intersection of multiple nodes."},

    /* Surface creation (with automatic registration for raycast) */
    {"sphere_surface", (PyCFunction)AleaTHORSystem_sphere_surface, METH_VARARGS,
     "sphere_surface(surface_id, cx, cy, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create sphere surface with both halfspace nodes registered for raycast."},
    {"cylinder_z_surface", (PyCFunction)AleaTHORSystem_cylinder_z_surface, METH_VARARGS,
     "cylinder_z_surface(surface_id, cx, cy, radius) -> (index, pos_node, neg_node)\n\n"
     "Create Z-cylinder surface with both halfspace nodes registered for raycast."},
    {"box_surface", (PyCFunction)AleaTHORSystem_box_surface, METH_VARARGS,
     "box_surface(surface_id, xmin, xmax, ymin, ymax, zmin, zmax) -> (index, pos_node, neg_node)\n\n"
     "Create box surface with both halfspace nodes registered for raycast."},
    {"plane_surface", (PyCFunction)AleaTHORSystem_plane_surface, METH_VARARGS,
     "plane_surface(surface_id, a, b, c, d) -> (index, pos_node, neg_node)\n\n"
     "Create plane surface with both halfspace nodes registered for raycast."},
    {"cylinder_x_surface", (PyCFunction)AleaTHORSystem_cylinder_x_surface, METH_VARARGS,
     "cylinder_x_surface(surface_id, cy, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create X-cylinder surface with both halfspace nodes registered for raycast."},
    {"cylinder_y_surface", (PyCFunction)AleaTHORSystem_cylinder_y_surface, METH_VARARGS,
     "cylinder_y_surface(surface_id, cx, cz, radius) -> (index, pos_node, neg_node)\n\n"
     "Create Y-cylinder surface with both halfspace nodes registered for raycast."},
    {"cone_z_surface", (PyCFunction)AleaTHORSystem_cone_z_surface, METH_VARARGS,
     "cone_z_surface(surface_id, cx, cy, cz, t_squared) -> (index, pos_node, neg_node)\n\n"
     "Create Z-cone surface with both halfspace nodes registered for raycast."},
    {"cone_x_surface", (PyCFunction)AleaTHORSystem_cone_x_surface, METH_VARARGS,
     "cone_x_surface(surface_id, cx, cy, cz, t_squared) -> (index, pos_node, neg_node)\n\n"
     "Create X-cone surface with both halfspace nodes registered for raycast."},
    {"cone_y_surface", (PyCFunction)AleaTHORSystem_cone_y_surface, METH_VARARGS,
     "cone_y_surface(surface_id, cx, cy, cz, t_squared) -> (index, pos_node, neg_node)\n\n"
     "Create Y-cone surface with both halfspace nodes registered for raycast."},
    {"torus_z_surface", (PyCFunction)AleaTHORSystem_torus_z_surface, METH_VARARGS,
     "torus_z_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create Z-torus surface with both halfspace nodes registered for raycast."},
    {"torus_x_surface", (PyCFunction)AleaTHORSystem_torus_x_surface, METH_VARARGS,
     "torus_x_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create X-torus surface with both halfspace nodes registered for raycast."},
    {"torus_y_surface", (PyCFunction)AleaTHORSystem_torus_y_surface, METH_VARARGS,
     "torus_y_surface(surface_id, cx, cy, cz, major_radius, minor_radius) -> (index, pos_node, neg_node)\n\n"
     "Create Y-torus surface with both halfspace nodes registered for raycast."},
    {"quadric_surface", (PyCFunction)AleaTHORSystem_quadric_surface, METH_VARARGS,
     "quadric_surface(surface_id, A, B, C, D, E, F, G, H, I, J) -> (index, pos_node, neg_node)\n\n"
     "Create general quadric surface (Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0)."},
    {"rcc_surface", (PyCFunction)AleaTHORSystem_rcc_surface, METH_VARARGS,
     "rcc_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, radius) -> (index, pos_node, neg_node)\n\n"
     "Create RCC (Right Circular Cylinder) macrobody surface."},
    {"box_general_surface", (PyCFunction)AleaTHORSystem_box_general_surface, METH_VARARGS,
     "box_general_surface(surface_id, corner_x, corner_y, corner_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z) -> (index, pos_node, neg_node)\n\n"
     "Create general box macrobody from corner and three edge vectors."},
    {"sph_surface", (PyCFunction)AleaTHORSystem_sph_surface, METH_VARARGS,
     "sph_surface(surface_id, cx, cy, cz, r) -> (index, pos_node, neg_node)\n\n"
     "Create SPH (sphere) macrobody surface."},
    {"trc_surface", (PyCFunction)AleaTHORSystem_trc_surface, METH_VARARGS,
     "trc_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, base_radius, top_radius) -> (index, pos_node, neg_node)\n\n"
     "Create TRC (Truncated Right Cone) macrobody surface."},
    {"ell_surface", (PyCFunction)AleaTHORSystem_ell_surface, METH_VARARGS,
     "ell_surface(surface_id, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, major_axis_len) -> (index, pos_node, neg_node)\n\n"
     "Create ELL (ellipsoid) macrobody from two foci and major axis length."},
    {"rec_surface", (PyCFunction)AleaTHORSystem_rec_surface, METH_VARARGS,
     "rec_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, axis1_x, axis1_y, axis1_z, axis2_x, axis2_y, axis2_z) -> (index, pos_node, neg_node)\n\n"
     "Create REC (Right Elliptical Cylinder) macrobody surface."},
    {"wed_surface", (PyCFunction)AleaTHORSystem_wed_surface, METH_VARARGS,
     "wed_surface(surface_id, vertex_x, vertex_y, vertex_z, v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z) -> (index, pos_node, neg_node)\n\n"
     "Create WED (wedge) macrobody from vertex and three edge vectors."},
    {"rhp_surface", (PyCFunction)AleaTHORSystem_rhp_surface, METH_VARARGS,
     "rhp_surface(surface_id, base_x, base_y, base_z, height_x, height_y, height_z, r1_x, r1_y, r1_z, r2_x, r2_y, r2_z, r3_x, r3_y, r3_z) -> (index, pos_node, neg_node)\n\n"
     "Create RHP (Right Hexagonal Prism) macrobody surface."},
    {"get_surface_nodes", (PyCFunction)AleaTHORSystem_get_surface_nodes, METH_VARARGS,
     "get_surface_nodes(surface_index) -> (pos_node, neg_node)\n\n"
     "Get halfspace nodes for a registered surface."},

    /* Cell registration */
    {"add_cell", (PyCFunction)AleaTHORSystem_add_cell,
     METH_VARARGS | METH_KEYWORDS,
     "add_cell(cell_id, root_node, material_id=0, density=0.0, universe_id=0) -> index"},

    /* Export */
    {"export_mcnp", (PyCFunction)AleaTHORSystem_export_mcnp,
     METH_VARARGS | METH_KEYWORDS,
     "export_mcnp(filename, deduplicate=True, universe_depth=-1, fill_depth=0)\n\n"
     "Export to MCNP format.\n\n"
     "Args:\n"
     "  filename: Output file path\n"
     "  deduplicate: Deduplicate surfaces (default True)\n"
     "  universe_depth: Universe filter (-1=all, 0=base only, N=N levels deep)\n"
     "  fill_depth: FILL expansion (0=none, N=N levels, -1=full flatten)"},
    {"export_openmc", (PyCFunction)AleaTHORSystem_export_openmc,
     METH_VARARGS | METH_KEYWORDS,
     "export_openmc(filename)\n\nExport to OpenMC XML format."},

    /* Merge */
    {"merge", (PyCFunction)AleaTHORSystem_merge, METH_VARARGS | METH_KEYWORDS,
     "merge(other, id_offset=0)\n\nMerge another system into this one."},

    /* Utilities */
    {"set_verbose", (PyCFunction)AleaTHORSystem_set_verbose, METH_VARARGS,
     "set_verbose(enabled)\n\nEnable/disable verbose output."},
    {"validate", (PyCFunction)AleaTHORSystem_validate, METH_NOARGS,
     "validate() -> int\n\nValidate system integrity."},
    {"print_summary", (PyCFunction)AleaTHORSystem_print_summary, METH_NOARGS,
     "print_summary()\n\nPrint system summary."},
    {"set_tolerance", (PyCFunction)AleaTHORSystem_set_tolerance,
     METH_VARARGS | METH_KEYWORDS,
     "set_tolerance(abs_tol=1e-6, rel_tol=1e-9, zero_thresh=1e-10)\n\nSet tolerances."},
    {"clone", (PyCFunction)AleaTHORSystem_clone, METH_NOARGS,
     "clone() -> System\n\nCreate deep copy of system."},
    {"reset", (PyCFunction)AleaTHORSystem_reset, METH_NOARGS,
     "reset()\n\nReset system to empty state."},

    /* Raycast */
    {"raycast", (PyCFunction)AleaTHORSystem_raycast, METH_VARARGS | METH_KEYWORDS,
     "raycast(ox, oy, oz, dx, dy, dz, t_max=0) -> list of segments\n\n"
     "Cast ray through geometry. Returns list of {t_enter, t_exit, cell_id, material_id, density}."},
    {"ray_first_cell", (PyCFunction)AleaTHORSystem_ray_first_cell, METH_VARARGS,
     "ray_first_cell(ox, oy, oz, dx, dy, dz, t_max=0) -> (cell_id, t) or None"},

    /* Slice curves API (for matplotlib) */
    {"get_slice_curves_z", (PyCFunction)AleaTHORSystem_get_slice_curves_z, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_z(z, x_min, x_max, y_min, y_max) -> dict\n\n"
     "Get analytical curves from surface-plane intersections (XY plane at z).\n\n"
     "Returns dict with:\n"
     "  'curves': list of curve dicts (type, surface_id, center/radius/etc)\n"
     "  'u_min', 'u_max', 'v_min', 'v_max': bounding box of curves"},
    {"get_slice_curves_y", (PyCFunction)AleaTHORSystem_get_slice_curves_y, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_y(y, x_min, x_max, z_min, z_max) -> dict\n\n"
     "Get analytical curves (XZ plane at y)."},
    {"get_slice_curves_x", (PyCFunction)AleaTHORSystem_get_slice_curves_x, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves_x(x, y_min, y_max, z_min, z_max) -> dict\n\n"
     "Get analytical curves (YZ plane at x)."},

    /* Grid cell queries (for matplotlib fills) */
    {"find_cells_grid_z", (PyCFunction)AleaTHORSystem_find_cells_grid_z, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_z(z, x_min, x_max, y_min, y_max, nx, ny, universe_depth=-1, detect_errors=False) -> dict\n\n"
     "Find cells on a 2D grid (XY plane at z).\n\n"
     "Args:\n"
     "  universe_depth: Universe depth (-1=innermost, 0=root, N=depth N)\n"
     "  detect_errors: Include error detection (overlaps, undefined)\n\n"
     "Returns dict with:\n"
     "  'cell_ids': list of cell IDs (nx*ny, -1 for void)\n"
     "  'material_ids': list of material IDs\n"
     "  'errors': list of error codes (if detect_errors=True): 0=ok, 1=overlap, 2=undefined\n"
     "  'nx', 'ny': grid dimensions"},
    {"find_cells_grid_y", (PyCFunction)AleaTHORSystem_find_cells_grid_y, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_y(y, x_min, x_max, z_min, z_max, nx, nz, universe_depth=-1, detect_errors=False) -> dict\n\n"
     "Find cells on XZ grid at y. Same optional args as find_cells_grid_z."},
    {"find_cells_grid_x", (PyCFunction)AleaTHORSystem_find_cells_grid_x, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid_x(x, y_min, y_max, z_min, z_max, ny, nz, universe_depth=-1, detect_errors=False) -> dict\n\n"
     "Find cells on YZ grid at x. Same optional args as find_cells_grid_z."},

    /* Arbitrary plane slice curves and grid */
    {"get_slice_curves", (PyCFunction)AleaTHORSystem_get_slice_curves, METH_VARARGS | METH_KEYWORDS,
     "get_slice_curves(origin, normal, up, u_min, u_max, v_min, v_max) -> dict\n\n"
     "Get analytical curves on arbitrary slice plane.\n\n"
     "Args:\n"
     "  origin: Point on the plane as (x, y, z) tuple\n"
     "  normal: Plane normal as (nx, ny, nz) tuple\n"
     "  up: Up vector hint as (ux, uy, uz) tuple\n"
     "  u_min, u_max: Horizontal bounds in plane coordinates\n"
     "  v_min, v_max: Vertical bounds in plane coordinates\n\n"
     "Returns dict with:\n"
     "  'curves': list of curve dicts\n"
     "  'u_min', 'u_max', 'v_min', 'v_max': bounds"},
    {"find_cells_grid", (PyCFunction)AleaTHORSystem_find_cells_grid, METH_VARARGS | METH_KEYWORDS,
     "find_cells_grid(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv, universe_depth=-1, detect_errors=False) -> dict\n\n"
     "Find cells on arbitrary slice plane grid.\n\n"
     "Args:\n"
     "  origin: Point on the plane as (x, y, z) tuple\n"
     "  normal: Plane normal as (nx, ny, nz) tuple\n"
     "  up: Up vector hint as (ux, uy, uz) tuple\n"
     "  u_min, u_max: Horizontal bounds in plane coordinates\n"
     "  v_min, v_max: Vertical bounds in plane coordinates\n"
     "  nu, nv: Grid resolution\n"
     "  universe_depth: Universe depth (-1=innermost, 0=root, N=depth N)\n"
     "  detect_errors: Include error detection (overlaps, undefined)\n\n"
     "Returns dict with:\n"
     "  'cell_ids': list of cell IDs (nu*nv, -1 for void)\n"
     "  'material_ids': list of material IDs\n"
     "  'errors': list of error codes (if detect_errors=True)\n"
     "  'nu', 'nv': grid dimensions"},

    /* Label positioning */
    {"find_label_positions", (PyCFunction)AleaTHORSystem_find_label_positions, METH_VARARGS | METH_KEYWORDS,
     "find_label_positions(ids, width, height, min_pixels=100) -> list\n\n"
     "Find optimal label positions for regions in a cell/material grid.\n\n"
     "Args:\n"
     "  ids: List of cell or material IDs from find_cells_grid_*()\n"
     "  width, height: Grid dimensions\n"
     "  min_pixels: Minimum region size to include (default 100)\n\n"
     "Returns list of dicts with:\n"
     "  'id': Cell/material ID\n"
     "  'px', 'py': Pixel coordinates for label\n"
     "  'pixel_count': Region size in pixels"},

    /* Cell filtering */
    {"get_cells_by_material", (PyCFunction)AleaTHORSystem_get_cells_by_material, METH_VARARGS,
     "get_cells_by_material(material_id) -> list of indices\n\n"
     "Get indices of all cells with the given material ID."},
    {"get_cells_by_universe", (PyCFunction)AleaTHORSystem_get_cells_by_universe, METH_VARARGS,
     "get_cells_by_universe(universe_id) -> list of indices\n\n"
     "Get indices of all cells in the given universe."},
    {"get_cells_filling_universe", (PyCFunction)AleaTHORSystem_get_cells_filling_universe, METH_VARARGS,
     "get_cells_filling_universe(universe_id) -> list of indices\n\n"
     "Get indices of all cells that FILL the given universe."},
    {"get_cells_in_bbox", (PyCFunction)AleaTHORSystem_get_cells_in_bbox, METH_VARARGS,
     "get_cells_in_bbox(x_min, x_max, y_min, y_max, z_min, z_max) -> list of indices\n\n"
     "Get indices of cells whose bounding boxes intersect the given region."},

    /* Extract operations */
    {"extract_universe", (PyCFunction)AleaTHORSystem_extract_universe, METH_VARARGS,
     "extract_universe(universe_id) -> System\n\n"
     "Extract a universe and all universes it references into a new system."},
    {"extract_region", (PyCFunction)AleaTHORSystem_extract_region, METH_VARARGS,
     "extract_region(x_min, x_max, y_min, y_max, z_min, z_max) -> System\n\n"
     "Extract cells in a bounding box region into a new system."},

    /* Material operations */
    {"create_mixture", (PyCFunction)AleaTHORSystem_create_mixture, METH_VARARGS | METH_KEYWORDS,
     "create_mixture(material_ids, fractions, new_id=0) -> int\n\n"
     "Create a mixture of materials.\n\n"
     "Args:\n"
     "  material_ids: List of material IDs to mix\n"
     "  fractions: List of fractions (normalized automatically)\n"
     "  new_id: ID for new mixture (0 for auto-assign)\n\n"
     "Returns: Assigned material ID"},

    /* Cell fill / ID operations */
    {"set_fill", (PyCFunction)AleaTHORSystem_set_fill, METH_VARARGS,
     "set_fill(cell_index, fill_universe, transform=0)\n\nSet fill universe for a cell."},
    {"get_cell_id", (PyCFunction)AleaTHORSystem_get_cell_id, METH_VARARGS,
     "get_cell_id(cell_index) -> int\n\nGet MCNP cell ID from cell index."},
    {"cells_in_universe", (PyCFunction)AleaTHORSystem_cells_in_universe, METH_VARARGS,
     "cells_in_universe(universe_id) -> list of indices\n\nGet cell indices in a universe."},
    {"find_cell_at", (PyCFunction)AleaTHORSystem_find_cell_at, METH_VARARGS,
     "find_cell_at(x, y, z) -> (cell_id, material_id) or None\n\nFind cell and get both cell ID and material."},

    /* Surface operations */
    {"surface_find", (PyCFunction)AleaTHORSystem_surface_find, METH_VARARGS,
     "surface_find(surface_id) -> int or None\n\nFind surface index by MCNP surface ID."},
    {"surface_node", (PyCFunction)AleaTHORSystem_surface_node, METH_VARARGS,
     "surface_node(surface_id, sense) -> node_id or None\n\nGet node by surface ID and sense (-1 or +1)."},

    /* CSG Node Inspection (individual) */
    {"node_operation", (PyCFunction)AleaTHORSystem_node_operation, METH_VARARGS,
     "node_operation(node_id) -> int\n\nGet operation type: 0=primitive, 1=union, 2=intersection, 3=difference, 4=complement."},
    {"node_left", (PyCFunction)AleaTHORSystem_node_left, METH_VARARGS,
     "node_left(node_id) -> node_id or None\n\nGet left child of a boolean node."},
    {"node_right", (PyCFunction)AleaTHORSystem_node_right, METH_VARARGS,
     "node_right(node_id) -> node_id or None\n\nGet right child of a boolean node."},
    {"node_primitive_type", (PyCFunction)AleaTHORSystem_node_primitive_type, METH_VARARGS,
     "node_primitive_type(node_id) -> int\n\nGet primitive type for a leaf node."},
    {"node_primitive_id", (PyCFunction)AleaTHORSystem_node_primitive_id, METH_VARARGS,
     "node_primitive_id(node_id) -> int or None\n\nGet primitive ID for a leaf node."},
    {"node_sense", (PyCFunction)AleaTHORSystem_node_sense, METH_VARARGS,
     "node_sense(node_id) -> int\n\nGet sense for a primitive node (+1 or -1)."},
    {"node_surface_id", (PyCFunction)AleaTHORSystem_node_surface_id, METH_VARARGS,
     "node_surface_id(node_id) -> int\n\nGet MCNP surface ID associated with a primitive node."},

    /* Renumbering */
    {"renumber_cells", (PyCFunction)AleaTHORSystem_renumber_cells, METH_VARARGS,
     "renumber_cells(start_id) -> int\n\nRenumber all cells starting from start_id."},
    {"renumber_surfaces", (PyCFunction)AleaTHORSystem_renumber_surfaces, METH_VARARGS,
     "renumber_surfaces(start_id) -> int\n\nRenumber all surfaces starting from start_id."},
    {"offset_cell_ids", (PyCFunction)AleaTHORSystem_offset_cell_ids, METH_VARARGS,
     "offset_cell_ids(offset)\n\nAdd offset to all cell IDs."},
    {"offset_surface_ids", (PyCFunction)AleaTHORSystem_offset_surface_ids, METH_VARARGS,
     "offset_surface_ids(offset)\n\nAdd offset to all surface IDs."},
    {"offset_material_ids", (PyCFunction)AleaTHORSystem_offset_material_ids, METH_VARARGS,
     "offset_material_ids(offset)\n\nAdd offset to all material IDs."},

    /* Split / Expand */
    {"split_union_cells", (PyCFunction)AleaTHORSystem_split_union_cells, METH_NOARGS,
     "split_union_cells() -> int\n\nSplit cells with top-level unions into multiple simpler cells.\n"
     "Returns number of new cells created."},
    {"expand_macrobodies", (PyCFunction)AleaTHORSystem_expand_macrobodies, METH_NOARGS,
     "expand_macrobodies() -> int\n\nExpand all macrobodies to primitive surfaces."},

    /* Volume estimation */
    {"compute_bounding_sphere", (PyCFunction)AleaTHORSystem_compute_bounding_sphere, METH_VARARGS,
     "compute_bounding_sphere(tol=1.0) -> (cx, cy, cz, radius)\n\n"
     "Compute a tight bounding sphere for the entire model."},
    {"estimate_cell_volumes", (PyCFunction)AleaTHORSystem_estimate_cell_volumes,
     METH_VARARGS | METH_KEYWORDS,
     "estimate_cell_volumes(ox, oy, oz, radius, n_rays) -> dict\n\n"
     "Estimate cell volumes using random ray tracing (Cauchy-Crofton method).\n"
     "Returns dict with 'volumes' and 'rel_errors' lists."},
    {"estimate_instance_volumes", (PyCFunction)AleaTHORSystem_estimate_instance_volumes, METH_VARARGS,
     "estimate_instance_volumes(n_rays) -> dict\n\n"
     "Estimate volumes per cell instance (spatial-index aware).\n"
     "Requires spatial index to be built. Returns dict with 'volumes' and 'rel_errors'."},
    {"remove_cells_by_volume", (PyCFunction)AleaTHORSystem_remove_cells_by_volume, METH_VARARGS,
     "remove_cells_by_volume(volumes, threshold) -> int\n\n"
     "Remove cells whose estimated volume is below threshold.\n"
     "Returns number of cells removed."},

    /* BBox tightening */
    {"tighten_cell_bbox", (PyCFunction)AleaTHORSystem_tighten_cell_bbox, METH_VARARGS,
     "tighten_cell_bbox(cell_index, tol=1.0) -> (xmin, xmax, ymin, ymax, zmin, zmax)\n\n"
     "Tighten a single cell's bounding box via interval arithmetic."},
    {"tighten_all_bboxes", (PyCFunction)AleaTHORSystem_tighten_all_bboxes, METH_VARARGS,
     "tighten_all_bboxes(tol=1.0) -> int\n\n"
     "Tighten all cell bounding boxes. Returns number of cells tightened."},

    /* Cell-aware raycast */
    {"raycast_cell_aware", (PyCFunction)AleaTHORSystem_raycast_cell_aware,
     METH_VARARGS | METH_KEYWORDS,
     "raycast_cell_aware(ox, oy, oz, dx, dy, dz, t_max=0) -> list of segments\n\n"
     "Cell-aware raycast using per-cell surface index. More efficient than global raycast."},

    /* Grid overlap check */
    {"check_grid_overlaps", (PyCFunction)AleaTHORSystem_check_grid_overlaps,
     METH_VARARGS | METH_KEYWORDS,
     "check_grid_overlaps(origin, normal, up, u_min, u_max, v_min, v_max, nu, nv, cell_ids, errors, universe_depth=-1) -> list\n\n"
     "Check grid for overlapping cells (comprehensive). Returns updated error list."},

    /* Surface label positions */
    {"find_surface_label_positions", (PyCFunction)AleaTHORSystem_find_surface_label_positions,
     METH_VARARGS | METH_KEYWORDS,
     "find_surface_label_positions(origin, normal, up, u_min, u_max, v_min, v_max, width, height, margin=20) -> list\n\n"
     "Find label positions for surfaces on a slice plane."},

    /* Config */
    {"get_config", (PyCFunction)AleaTHORSystem_get_config, METH_NOARGS,
     "get_config() -> dict\n\nGet current system configuration."},
    {"set_config", (PyCFunction)AleaTHORSystem_set_config, METH_VARARGS,
     "set_config(dict)\n\nUpdate system configuration from dict. Only provided keys are changed."},
    {"set_log_level", (PyCFunction)AleaTHORSystem_set_log_level, METH_VARARGS,
     "set_log_level(level)\n\nSet log level: 0=none, 1=error, 2=warn, 3=info, 4=debug, 5=trace."},

    /* CSG Simplification */
    {"flatten_all_cells", (PyCFunction)AleaTHORSystem_flatten_all_cells, METH_NOARGS,
     "flatten_all_cells() -> dict\n\nFull CSG simplification pass on all cells.\n"
     "Returns dict with simplification statistics."},

    /* Numerical BBox tightening */
    {"tighten_cell_bbox_numerical", (PyCFunction)AleaTHORSystem_tighten_cell_bbox_numerical, METH_VARARGS,
     "tighten_cell_bbox_numerical(cell_index) -> None\n\n"
     "Tighten a cell's bounding box using numerical sampling (fallback for complex cells)."},

    /* Primitive data */
    {"node_primitive_data", (PyCFunction)AleaTHORSystem_node_primitive_data, METH_VARARGS,
     "node_primitive_data(node_id) -> dict\n\n"
     "Get full primitive geometry data for a leaf node.\n"
     "Returns dict with 'type' and type-specific fields (center, radius, coefficients, etc.)."},

    /* Mesh module */
    {"mesh_export", (PyCFunction)AleaTHORSystem_mesh_export, METH_VARARGS | METH_KEYWORDS,
     "mesh_export(filename, nx=10, ny=10, nz=10, ...) -> None\n\n"
     "Export geometry as structured mesh (Gmsh or VTK).\n"
     "Optional kwargs: x_min, x_max, y_min, y_max, z_min, z_max (auto if 0),\n"
     "format ('gmsh'/'vtk'), void_material_id, auto_pad."},
    {"mesh_sample", (PyCFunction)AleaTHORSystem_mesh_sample, METH_VARARGS | METH_KEYWORDS,
     "mesh_sample(nx=10, ny=10, nz=10, ...) -> dict\n\n"
     "Sample geometry on structured mesh. Returns dict with:\n"
     "  'material_ids', 'cell_ids': flat lists (nx*ny*nz, Z-major)\n"
     "  'x_nodes', 'y_nodes', 'z_nodes': node positions (n+1 each)\n"
     "  'nx', 'ny', 'nz': grid dimensions."},

    {NULL}
};


/* ============================================================================
 * AleaTHORSystem Type Definition
 * ============================================================================ */

static PyTypeObject AleaTHORSystemType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.System",
    .tp_doc = PyDoc_STR("AleaTHOR-CSG geometry system.\n\n"
                        "Create an empty system with System(), or use load_mcnp() to load a file."),
    .tp_basicsize = sizeof(AleaTHORSystemObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = AleaTHORSystem_new,
    .tp_init = (initproc)AleaTHORSystem_init,
    .tp_dealloc = (destructor)AleaTHORSystem_dealloc,
    .tp_methods = AleaTHORSystem_methods,
    .tp_getset = AleaTHORSystem_getsetters,
};

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

/* ============================================================================
 * Module-level Functions
 * ============================================================================ */

static PyObject* mod_load_mcnp(PyObject* self, PyObject* args) {
    (void)self;
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

    alea_system_t* sys;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    sys = alea_load_mcnp(filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (sys) alea_destroy(sys);
        return NULL;
    }

    if (!sys) {
        PyErr_Format(PyExc_IOError, "Failed to load %s: %s", filename, alea_error());
        return NULL;
    }

    AleaTHORSystemObject* obj = (AleaTHORSystemObject*)AleaTHORSystemType.tp_alloc(&AleaTHORSystemType, 0);
    if (!obj) {
        alea_destroy(sys);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 1;
    return (PyObject*)obj;
}

static PyObject* mod_load_openmc(PyObject* self, PyObject* args) {
    (void)self;
    const char* filename;
    if (!PyArg_ParseTuple(args, "s", &filename)) return NULL;

    alea_system_t* sys;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    sys = alea_load_openmc(filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (sys) alea_destroy(sys);
        return NULL;
    }

    if (!sys) {
        PyErr_Format(PyExc_IOError, "Failed to load %s: %s", filename, alea_error());
        return NULL;
    }

    AleaTHORSystemObject* obj = (AleaTHORSystemObject*)AleaTHORSystemType.tp_alloc(&AleaTHORSystemType, 0);
    if (!obj) {
        alea_destroy(sys);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 1;
    return (PyObject*)obj;
}

static PyObject* mod_load_mcnp_string(PyObject* self, PyObject* args) {
    (void)self;
    const char* input;
    Py_ssize_t length;

    if (!PyArg_ParseTuple(args, "s#", &input, &length)) return NULL;

    alea_system_t* sys;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    sys = alea_load_mcnp_string(input, (size_t)length);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (sys) alea_destroy(sys);
        return NULL;
    }

    if (!sys) {
        PyErr_Format(PyExc_ValueError, "Failed to parse MCNP input: %s", alea_error());
        return NULL;
    }

    AleaTHORSystemObject* obj = (AleaTHORSystemObject*)AleaTHORSystemType.tp_alloc(&AleaTHORSystemType, 0);
    if (!obj) {
        alea_destroy(sys);
        return NULL;
    }

    obj->sys = sys;
    obj->owns_sys = 1;
    return (PyObject*)obj;
}

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

static PyObject* mod_version(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    return PyUnicode_FromString(alea_version());
}

static PyObject* mod_get_error(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    const char* err = alea_error();
    if (err && *err) return PyUnicode_FromString(err);
    Py_RETURN_NONE;
}

static PyObject* mod_clear_error(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    alea_error_clear();
    Py_RETURN_NONE;
}

/* ============================================================================
 * Logging Integration with Python's logging module
 * ============================================================================ */

/* Cached Python logging module and logger */
static PyObject* g_logging_module = NULL;
static PyObject* g_logger = NULL;

/* Map C log levels to Python logging method names */
static const char* level_to_method_name(alea_log_level_t level) {
    switch (level) {
        case ALEA_LOG_LEVEL_ERROR: return "error";
        case ALEA_LOG_LEVEL_WARN:  return "warning";
        case ALEA_LOG_LEVEL_INFO:  return "info";
        case ALEA_LOG_LEVEL_DEBUG: return "debug";
        case ALEA_LOG_LEVEL_TRACE: return "debug";  /* Python has no TRACE, use DEBUG */
        default:            return "debug";
    }
}

/* Callback that routes C library logs to Python logging module */
static void python_log_callback(alea_log_level_t level, const char* file,
                                 int line, const char* message, void* user_data) {
    (void)user_data;

    /* Must hold GIL to call Python */
    PyGILState_STATE gstate = PyGILState_Ensure();

    if (g_logger) {
        /* Build log message with file:line for debug levels */
        PyObject* log_msg;
        if (file && (level >= ALEA_LOG_LEVEL_DEBUG)) {
            log_msg = PyUnicode_FromFormat("[%s:%d] %s", file, line, message);
        } else {
            log_msg = PyUnicode_FromString(message);
        }

        if (log_msg) {
            const char* method = level_to_method_name(level);
            PyObject* result = PyObject_CallMethod(g_logger, method, "O", log_msg);
            Py_XDECREF(result);
            Py_DECREF(log_msg);
        }

        /* Clear any Python exception (don't let logging errors propagate) */
        if (PyErr_Occurred()) {
            PyErr_Clear();
        }
    }

    PyGILState_Release(gstate);
}

/* Initialize logging integration */
static int init_logging(void) {
    /* Import logging module */
    g_logging_module = PyImport_ImportModule("logging");
    if (!g_logging_module) {
        PyErr_Clear();
        return -1;
    }

    /* Get a logger named "aleathor" */
    g_logger = PyObject_CallMethod(g_logging_module, "getLogger", "s", "aleathor");
    if (!g_logger) {
        PyErr_Clear();
        Py_CLEAR(g_logging_module);
        return -1;
    }

    /* Set C library to emit ALL messages - Python logging handles filtering */
    alea_log_set_level(ALEA_LOG_LEVEL_TRACE);

    /* Set the C library callback to route to Python */
    alea_log_set_callback(python_log_callback, NULL);

    return 0;
}

/* Cleanup logging on module unload */
static void cleanup_logging(void) {
    alea_log_set_callback(NULL, NULL);
    Py_CLEAR(g_logger);
    Py_CLEAR(g_logging_module);
}

/* Map our log levels to Python logging levels */
static int to_python_log_level(int level) {
    switch (level) {
        case ALEA_LOG_LEVEL_NONE:  return 100;  /* Higher than CRITICAL, effectively disabled */
        case ALEA_LOG_LEVEL_ERROR: return 40;   /* logging.ERROR */
        case ALEA_LOG_LEVEL_WARN:  return 30;   /* logging.WARNING */
        case ALEA_LOG_LEVEL_INFO:  return 20;   /* logging.INFO */
        case ALEA_LOG_LEVEL_DEBUG: return 10;   /* logging.DEBUG */
        case ALEA_LOG_LEVEL_TRACE: return 5;    /* Below DEBUG */
        default:            return 30;   /* Default to WARNING */
    }
}

/* Map Python logging level back to our levels */
static int from_python_log_level(int py_level) {
    if (py_level >= 100) return ALEA_LOG_LEVEL_NONE;
    if (py_level >= 40) return ALEA_LOG_LEVEL_ERROR;
    if (py_level >= 30) return ALEA_LOG_LEVEL_WARN;
    if (py_level >= 20) return ALEA_LOG_LEVEL_INFO;
    if (py_level >= 10) return ALEA_LOG_LEVEL_DEBUG;
    return ALEA_LOG_LEVEL_TRACE;
}

/* Python API: set_log_level(level) - sets the Python logger level */
static PyObject* mod_set_log_level(PyObject* self, PyObject* args) {
    (void)self;
    int level;
    if (!PyArg_ParseTuple(args, "i", &level)) return NULL;

    if (level < ALEA_LOG_LEVEL_NONE || level > ALEA_LOG_LEVEL_TRACE) {
        PyErr_SetString(PyExc_ValueError, "Invalid log level (must be 0-5)");
        return NULL;
    }

    if (g_logger) {
        int py_level = to_python_log_level(level);
        PyObject* result = PyObject_CallMethod(g_logger, "setLevel", "i", py_level);
        Py_XDECREF(result);
        if (PyErr_Occurred()) {
            return NULL;
        }
    }

    Py_RETURN_NONE;
}

/* Python API: get_log_level() -> int */
static PyObject* mod_get_log_level(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;

    if (g_logger) {
        PyObject* effective_level = PyObject_CallMethod(g_logger, "getEffectiveLevel", NULL);
        if (effective_level && PyLong_Check(effective_level)) {
            int py_level = (int)PyLong_AsLong(effective_level);
            Py_DECREF(effective_level);
            return PyLong_FromLong(from_python_log_level(py_level));
        }
        Py_XDECREF(effective_level);
        PyErr_Clear();
    }

    return PyLong_FromLong(ALEA_LOG_LEVEL_WARN);  /* Default */
}

/* Python API: disable_logging() - restore default stderr output */
static PyObject* mod_disable_logging(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    alea_log_set_callback(NULL, NULL);
    Py_RETURN_NONE;
}

/* Python API: enable_logging() - route to Python logging module */
static PyObject* mod_enable_logging(PyObject* self, PyObject* Py_UNUSED(ignored)) {
    (void)self;
    if (!g_logger) {
        if (init_logging() < 0) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize Python logging integration");
            return NULL;
        }
    }
    alea_log_set_callback(python_log_callback, NULL);
    Py_RETURN_NONE;
}

/* ============================================================================
 * Module Definition
 * ============================================================================ */

static PyMethodDef mod_methods[] = {
    {"load_mcnp", mod_load_mcnp, METH_VARARGS,
     "load_mcnp(filename) -> System\n\nLoad MCNP input file."},
    {"load_mcnp_string", mod_load_mcnp_string, METH_VARARGS,
     "load_mcnp_string(input) -> System\n\nLoad MCNP from string."},
    {"load_openmc", mod_load_openmc, METH_VARARGS,
     "load_openmc(filename) -> System\n\nLoad OpenMC XML geometry file."},
    {"generate_void", (PyCFunction)mod_generate_void, METH_VARARGS | METH_KEYWORDS,
     "generate_void(system, bounds=None, max_depth=8, min_size=0.1, samples_per_node=27) -> VoidResult\n\n"
     "Generate void regions using octree algorithm."},
    {"version", mod_version, METH_NOARGS,
     "version() -> str\n\nGet library version string."},
    {"get_error", mod_get_error, METH_NOARGS,
     "get_error() -> str or None\n\nGet last error message."},
    {"clear_error", mod_clear_error, METH_NOARGS,
     "clear_error()\n\nClear error state."},
    /* Logging */
    {"set_log_level", mod_set_log_level, METH_VARARGS,
     "set_log_level(level)\n\n"
     "Set log level for the 'aleathor' logger.\n\n"
     "Levels: LOG_NONE=0, LOG_ERROR=1, LOG_WARN=2, LOG_INFO=3, LOG_DEBUG=4, LOG_TRACE=5\n\n"
     "Example:\n"
     "    import aleathor as ath\n"
     "    ath.set_log_level(ath.LOG_DEBUG)  # See debug messages"},
    {"get_log_level", mod_get_log_level, METH_NOARGS,
     "get_log_level() -> int\n\nGet current log level of the 'aleathor' logger."},
    {"enable_logging", mod_enable_logging, METH_NOARGS,
     "enable_logging()\n\n"
     "Route C library logs to Python's logging module (logger: 'aleathor').\n"
     "This is enabled by default on module import."},
    {"disable_logging", mod_disable_logging, METH_NOARGS,
     "disable_logging()\n\n"
     "Disable Python logging integration, restore default stderr output."},
    {NULL, NULL, 0, NULL}
};

/* Module cleanup function */
static void mod_module_free(void* module) {
    (void)module;
    cleanup_logging();
}

static struct PyModuleDef mod_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_alea",
    .m_doc = PyDoc_STR("AleaTHOR-CSG: Constructive Solid Geometry library for nuclear simulations.\n\n"
                       "This module provides Python bindings for the AleaTHOR-CSG library,\n"
                       "which parses MCNP input files and provides CSG tree operations.\n\n"
                       "Logging is automatically routed to Python's logging module (logger: 'aleathor')."),
    .m_size = 0,  /* Use 0 instead of -1 to enable m_free */
    .m_methods = mod_methods,
    .m_free = mod_module_free,
};

PyMODINIT_FUNC PyInit__alea(void) {
    /* Initialize types */
    if (PyType_Ready(&AleaTHORSystemType) < 0) return NULL;
    if (PyType_Ready(&AleaTHORVoidResultType) < 0) return NULL;

    /* Create module */
    PyObject* m = PyModule_Create(&mod_module);
    if (!m) return NULL;

    /* Add types */
    Py_INCREF(&AleaTHORSystemType);
    if (PyModule_AddObject(m, "System", (PyObject*)&AleaTHORSystemType) < 0) {
        Py_DECREF(&AleaTHORSystemType);
        Py_DECREF(m);
        return NULL;
    }

    Py_INCREF(&AleaTHORVoidResultType);
    if (PyModule_AddObject(m, "VoidResult", (PyObject*)&AleaTHORVoidResultType) < 0) {
        Py_DECREF(&AleaTHORVoidResultType);
        Py_DECREF(m);
        return NULL;
    }

    /* Add constants */
    PyModule_AddIntConstant(m, "NODE_INVALID", ALEA_NODE_ID_INVALID);
    PyModule_AddIntConstant(m, "MATERIAL_NONE", ALEA_MATERIAL_NONE);

    /* Add log level constants */
    PyModule_AddIntConstant(m, "LOG_NONE", ALEA_LOG_LEVEL_NONE);
    PyModule_AddIntConstant(m, "LOG_ERROR", ALEA_LOG_LEVEL_ERROR);
    PyModule_AddIntConstant(m, "LOG_WARN", ALEA_LOG_LEVEL_WARN);
    PyModule_AddIntConstant(m, "LOG_INFO", ALEA_LOG_LEVEL_INFO);
    PyModule_AddIntConstant(m, "LOG_DEBUG", ALEA_LOG_LEVEL_DEBUG);
    PyModule_AddIntConstant(m, "LOG_TRACE", ALEA_LOG_LEVEL_TRACE);

    /* Add CSG operation types */
    PyModule_AddIntConstant(m, "ALEA_OP_PRIMITIVE", ALEA_OP_PRIMITIVE);
    PyModule_AddIntConstant(m, "ALEA_OP_UNION", ALEA_OP_UNION);
    PyModule_AddIntConstant(m, "ALEA_OP_INTERSECTION", ALEA_OP_INTERSECTION);
    PyModule_AddIntConstant(m, "ALEA_OP_DIFFERENCE", ALEA_OP_DIFFERENCE);
    PyModule_AddIntConstant(m, "ALEA_OP_COMPLEMENT", ALEA_OP_COMPLEMENT);

    /* Add primitive types */
    PyModule_AddIntConstant(m, "PRIMITIVE_PLANE", ALEA_PRIMITIVE_PLANE);
    PyModule_AddIntConstant(m, "PRIMITIVE_SPHERE", ALEA_PRIMITIVE_SPHERE);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_X", ALEA_PRIMITIVE_CYLINDER_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_Y", ALEA_PRIMITIVE_CYLINDER_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_CYLINDER_Z", ALEA_PRIMITIVE_CYLINDER_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_X", ALEA_PRIMITIVE_CONE_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_Y", ALEA_PRIMITIVE_CONE_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_CONE_Z", ALEA_PRIMITIVE_CONE_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_RPP", ALEA_PRIMITIVE_RPP);
    PyModule_AddIntConstant(m, "PRIMITIVE_QUADRIC", ALEA_PRIMITIVE_QUADRIC);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_X", ALEA_PRIMITIVE_TORUS_X);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_Y", ALEA_PRIMITIVE_TORUS_Y);
    PyModule_AddIntConstant(m, "PRIMITIVE_TORUS_Z", ALEA_PRIMITIVE_TORUS_Z);
    PyModule_AddIntConstant(m, "PRIMITIVE_RCC", ALEA_PRIMITIVE_RCC);
    PyModule_AddIntConstant(m, "PRIMITIVE_BOX", ALEA_PRIMITIVE_BOX);
    PyModule_AddIntConstant(m, "PRIMITIVE_SPH", ALEA_PRIMITIVE_SPH);
    PyModule_AddIntConstant(m, "PRIMITIVE_TRC", ALEA_PRIMITIVE_TRC);
    PyModule_AddIntConstant(m, "PRIMITIVE_ELL", ALEA_PRIMITIVE_ELL);
    PyModule_AddIntConstant(m, "PRIMITIVE_REC", ALEA_PRIMITIVE_REC);
    PyModule_AddIntConstant(m, "PRIMITIVE_WED", ALEA_PRIMITIVE_WED);
    PyModule_AddIntConstant(m, "PRIMITIVE_RHP", ALEA_PRIMITIVE_RHP);
    PyModule_AddIntConstant(m, "PRIMITIVE_ARB", ALEA_PRIMITIVE_ARB);

    /* Add version info */
    PyModule_AddIntConstant(m, "VERSION_MAJOR", ALEA_VERSION_MAJOR);
    PyModule_AddIntConstant(m, "VERSION_MINOR", ALEA_VERSION_MINOR);
    PyModule_AddIntConstant(m, "VERSION_PATCH", ALEA_VERSION_PATCH);

    /* Initialize logging integration (route C logs to Python logging module) */
    init_logging();

    return m;
}
