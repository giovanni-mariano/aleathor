// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Slice curves (Z/Y/X/arbitrary), grid cell queries,
 *           label positions, grid overlap check, surface label positions.
 */

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
