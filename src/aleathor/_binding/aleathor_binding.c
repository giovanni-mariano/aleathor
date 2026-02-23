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
 *
 * Implementation is split across multiple files that are #included below
 * to keep a single compilation unit (required by the CPython extension API).
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

/* ============================================================================
 * Shared Helpers
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

/* ============================================================================
 * Implementation Files
 * ============================================================================ */

#include "_bind_core.c"
#include "_bind_geometry.c"
#include "_bind_io.c"
#include "_bind_util.c"
#include "_bind_raycast.c"
#include "_bind_slice.c"
#include "_bind_inspect.c"
#include "_bind_mesh.c"

/* ============================================================================
 * AleaTHORSystem Getters/Setters Table
 * ============================================================================ */

static PyGetSetDef AleaTHORSystem_getsetters[] = {
    {"cell_count", (getter)AleaTHORSystem_get_cell_count, NULL, "Number of cells", NULL},
    {"surface_count", (getter)AleaTHORSystem_get_surface_count, NULL, "Number of surfaces", NULL},
    {"universe_count", (getter)AleaTHORSystem_get_universe_count, NULL, "Number of universes", NULL},
    {"spatial_index_instance_count", (getter)AleaTHORSystem_get_spatial_index_instance_count, NULL, "Number of cell instances in spatial index", NULL},
    {NULL}
};

/* ============================================================================
 * AleaTHORSystem Method Table
 * ============================================================================ */

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
 * VoidResult Type + Module Init (included files)
 * ============================================================================ */

#include "_bind_void.c"
#include "_bind_module.c"
