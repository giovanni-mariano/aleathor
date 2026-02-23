// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Mesh export, mesh sample.
 */

/* ============================================================================
 * Mesh Module
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
