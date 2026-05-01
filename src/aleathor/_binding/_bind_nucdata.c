// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of aleathor_binding.c — do NOT compile separately.
 * It is #included from aleathor_binding.c to keep a single translation unit.
 *
 * Contents: Nuclear data Python types (XsDir, Nuclide, NucMaterial, Multigroup)
 *           and module-level functions for the nucdata API.
 */

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */

static PyTypeObject AleaTHORXsDirType;
static PyTypeObject AleaTHORNuclideType;
static PyTypeObject AleaTHORNucMaterialType;
static PyTypeObject AleaTHORMultigroupType;

/* ============================================================================
 * XsDir Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_xsdir_t* xsdir;
} AleaTHORXsDirObject;

static void AleaTHORXsDir_dealloc(AleaTHORXsDirObject* self) {
    if (self->xsdir) {
        alea_nuc_xsdir_free(self->xsdir);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORXsDir_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    AleaTHORXsDirObject* self = (AleaTHORXsDirObject*)type->tp_alloc(type, 0);
    if (self) self->xsdir = NULL;
    return (PyObject*)self;
}

static int AleaTHORXsDir_init(AleaTHORXsDirObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"path", "directory", NULL};
    const char* path = NULL;
    int directory = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|p", kwlist, &path, &directory))
        return -1;

    if (directory) {
        self->xsdir = alea_nuc_xsdir_load_dir(path);
    } else {
        self->xsdir = alea_nuc_xsdir_load(path);
    }

    if (!self->xsdir) {
        PyErr_Format(PyExc_IOError, "Failed to load xsdir from '%s'", path);
        return -1;
    }

    return 0;
}

static PyObject* AleaTHORXsDir_find(AleaTHORXsDirObject* self, PyObject* args) {
    const char* zaid;
    if (!PyArg_ParseTuple(args, "s", &zaid)) return NULL;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }

    const alea_nuc_xsdir_entry_t* entry = alea_nuc_xsdir_find(self->xsdir, zaid);
    if (!entry) Py_RETURN_NONE;

    return Py_BuildValue("{s:s, s:d, s:s, s:i, s:i, s:d}",
        "zaid", entry->zaid,
        "awr", entry->awr,
        "filename", entry->filename,
        "file_type", entry->file_type,
        "address", entry->address,
        "temperature", entry->temperature);
}

static PyObject* AleaTHORXsDir_get_count(AleaTHORXsDirObject* self, void* closure) {
    (void)closure;
    if (!self->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return NULL;
    }
    return PyLong_FromSize_t(alea_nuc_xsdir_count(self->xsdir));
}

static PyGetSetDef AleaTHORXsDir_getsetters[] = {
    {"count", (getter)AleaTHORXsDir_get_count, NULL, "Number of entries", NULL},
    {NULL}
};

static PyMethodDef AleaTHORXsDir_methods[] = {
    {"find", (PyCFunction)AleaTHORXsDir_find, METH_VARARGS,
     "find(zaid) -> dict or None\n\nFind xsdir entry by ZAID string (e.g. '92235.80c')."},
    {NULL}
};

static PyTypeObject AleaTHORXsDirType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.XsDir",
    .tp_doc = PyDoc_STR("Nuclear data cross-section directory (xsdir/xsdata)."),
    .tp_basicsize = sizeof(AleaTHORXsDirObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = AleaTHORXsDir_new,
    .tp_init = (initproc)AleaTHORXsDir_init,
    .tp_dealloc = (destructor)AleaTHORXsDir_dealloc,
    .tp_methods = AleaTHORXsDir_methods,
    .tp_getset = AleaTHORXsDir_getsetters,
};

/* ============================================================================
 * Nuclide Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_nuclide_t* nuc;
    int owned;                    /* 1 = we own it (must free), 0 = borrowed */
    PyObject* xsdir_ref;         /* Keep xsdir alive if borrowed */
} AleaTHORNuclideObject;

static void AleaTHORNuclide_dealloc(AleaTHORNuclideObject* self) {
    if (self->nuc && self->owned) {
        alea_nuc_nuclide_free(self->nuc);
    }
    Py_XDECREF(self->xsdir_ref);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORNuclide_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    AleaTHORNuclideObject* self = (AleaTHORNuclideObject*)type->tp_alloc(type, 0);
    if (self) {
        self->nuc = NULL;
        self->owned = 0;
        self->xsdir_ref = NULL;
    }
    return (PyObject*)self;
}

static int AleaTHORNuclide_init(AleaTHORNuclideObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"xsdir", "zaid", "cached", NULL};
    AleaTHORXsDirObject* xsdir_obj;
    const char* zaid;
    int cached = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!s|p", kwlist,
                                     &AleaTHORXsDirType, &xsdir_obj, &zaid, &cached))
        return -1;

    if (!xsdir_obj->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "XsDir not initialized");
        return -1;
    }

    if (cached) {
        /* Borrowed from xsdir cache */
        self->nuc = alea_nuc_xsdir_get_nuclide(xsdir_obj->xsdir, zaid);
        if (!self->nuc) {
            PyErr_Format(PyExc_ValueError, "Failed to load nuclide '%s'", zaid);
            return -1;
        }
        self->owned = 0;
        self->xsdir_ref = (PyObject*)xsdir_obj;
        Py_INCREF(self->xsdir_ref);
    } else {
        /* Owned copy */
        self->nuc = alea_nuc_load_nuclide(xsdir_obj->xsdir, zaid);
        if (!self->nuc) {
            PyErr_Format(PyExc_ValueError, "Failed to load nuclide '%s'", zaid);
            return -1;
        }
        self->owned = 1;
        self->xsdir_ref = NULL;
    }

    return 0;
}

/* Cross-section lookups */
static PyObject* AleaTHORNuclide_xs_total(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_total(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_xs_absorption(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_absorption(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_xs_elastic(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_elastic(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_xs_reaction(AleaTHORNuclideObject* self, PyObject* args) {
    int mt;
    double energy;
    if (!PyArg_ParseTuple(args, "id", &mt, &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_reaction(self->nuc, mt, energy));
}

static PyObject* AleaTHORNuclide_xs_heating(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_xs_heating(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_heating_per_collision(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_heating_per_collision(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_nu_bar(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->nuc->fission) {
        PyErr_SetString(PyExc_ValueError, "Nuclide is not fissile");
        return NULL;
    }
    return PyFloat_FromDouble(alea_nuc_nu_bar(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_doppler_broaden(AleaTHORNuclideObject* self, PyObject* args) {
    double kT_target;
    if (!PyArg_ParseTuple(args, "d", &kT_target)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    if (!self->owned) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot Doppler-broaden a cached (borrowed) nuclide. Load with cached=False.");
        return NULL;
    }
    alea_error_t err = alea_nuc_doppler_broaden(self->nuc, kT_target);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError, "Doppler broadening failed: %s", alea_error_string(err));
        return NULL;
    }
    Py_RETURN_NONE;
}

/* Photon cross sections */
static PyObject* AleaTHORNuclide_photon_xs_incoherent(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_incoherent(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_photon_xs_coherent(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_coherent(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_photon_xs_photoelectric(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_photoelectric(self->nuc, energy));
}

static PyObject* AleaTHORNuclide_photon_xs_pair(AleaTHORNuclideObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_photon_xs_pair(self->nuc, energy));
}

/* URR factors */
static PyObject* AleaTHORNuclide_urr_factors(AleaTHORNuclideObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    double factors[5];
    int applies = alea_nuc_urr_factors(self->nuc, energy, xi, factors);
    if (!applies) Py_RETURN_NONE;

    return Py_BuildValue("{s:d, s:d, s:d, s:d, s:d}",
        "total", factors[0], "elastic", factors[1],
        "fission", factors[2], "capture", factors[3],
        "heating", factors[4]);
}

/* Energy grid access */
static PyObject* AleaTHORNuclide_get_energy_grid(AleaTHORNuclideObject* self, PyObject* args) {
    (void)args;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    PyObject* list = PyList_New(self->nuc->n_energies);
    if (!list) return NULL;
    for (int i = 0; i < self->nuc->n_energies; i++) {
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(self->nuc->energy[i]));
    }
    return list;
}

/* Reactions list */
static PyObject* AleaTHORNuclide_get_reactions(AleaTHORNuclideObject* self, PyObject* args) {
    (void)args;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    PyObject* list = PyList_New(self->nuc->n_reactions);
    if (!list) return NULL;
    for (int i = 0; i < self->nuc->n_reactions; i++) {
        const alea_nuc_reaction_t* rxn = &self->nuc->reactions[i];
        PyObject* d = Py_BuildValue("{s:i, s:d, s:i}",
            "mt", rxn->mt,
            "q_value", rxn->q_value,
            "ty", rxn->ty);
        if (!d) { Py_DECREF(list); return NULL; }
        PyList_SET_ITEM(list, i, d);
    }
    return list;
}

/* Reaction yield and classification */
static PyObject* AleaTHORNuclide_reaction_yield(AleaTHORNuclideObject* self, PyObject* args) {
    int mt;
    double energy;
    if (!PyArg_ParseTuple(args, "id", &mt, &energy)) return NULL;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_reaction_yield(self->nuc, mt, energy));
}

/* Properties */
static PyObject* AleaTHORNuclide_get_zaid(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyUnicode_FromString(self->nuc->zaid);
}

static PyObject* AleaTHORNuclide_get_Z(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->Z);
}

static PyObject* AleaTHORNuclide_get_A(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->A);
}

static PyObject* AleaTHORNuclide_get_awr(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(self->nuc->awr);
}

static PyObject* AleaTHORNuclide_get_temperature(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyFloat_FromDouble(self->nuc->temperature);
}

static PyObject* AleaTHORNuclide_get_is_fissile(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->fission != NULL);
}

static PyObject* AleaTHORNuclide_get_has_urr(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->urr != NULL);
}

static PyObject* AleaTHORNuclide_get_has_photon(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyBool_FromLong(self->nuc->photon != NULL);
}

static PyObject* AleaTHORNuclide_get_n_reactions(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->n_reactions);
}

static PyObject* AleaTHORNuclide_get_n_energies(AleaTHORNuclideObject* self, void* closure) {
    (void)closure;
    if (!self->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }
    return PyLong_FromLong(self->nuc->n_energies);
}

static PyGetSetDef AleaTHORNuclide_getsetters[] = {
    {"zaid", (getter)AleaTHORNuclide_get_zaid, NULL, "ZAID string (e.g. '92235.80c')", NULL},
    {"Z", (getter)AleaTHORNuclide_get_Z, NULL, "Atomic number", NULL},
    {"A", (getter)AleaTHORNuclide_get_A, NULL, "Mass number", NULL},
    {"awr", (getter)AleaTHORNuclide_get_awr, NULL, "Atomic weight ratio", NULL},
    {"temperature", (getter)AleaTHORNuclide_get_temperature, NULL, "Temperature kT in MeV", NULL},
    {"is_fissile", (getter)AleaTHORNuclide_get_is_fissile, NULL, "True if nuclide has fission data", NULL},
    {"has_urr", (getter)AleaTHORNuclide_get_has_urr, NULL, "True if nuclide has URR probability tables", NULL},
    {"has_photon", (getter)AleaTHORNuclide_get_has_photon, NULL, "True if nuclide has photon data", NULL},
    {"n_reactions", (getter)AleaTHORNuclide_get_n_reactions, NULL, "Number of non-elastic reactions", NULL},
    {"n_energies", (getter)AleaTHORNuclide_get_n_energies, NULL, "Number of energy grid points", NULL},
    {NULL}
};

static PyMethodDef AleaTHORNuclide_methods[] = {
    {"xs_total", (PyCFunction)AleaTHORNuclide_xs_total, METH_VARARGS,
     "xs_total(energy) -> float\n\nTotal cross section in barns at energy (MeV)."},
    {"xs_absorption", (PyCFunction)AleaTHORNuclide_xs_absorption, METH_VARARGS,
     "xs_absorption(energy) -> float\n\nAbsorption cross section in barns."},
    {"xs_elastic", (PyCFunction)AleaTHORNuclide_xs_elastic, METH_VARARGS,
     "xs_elastic(energy) -> float\n\nElastic scattering cross section in barns."},
    {"xs_reaction", (PyCFunction)AleaTHORNuclide_xs_reaction, METH_VARARGS,
     "xs_reaction(mt, energy) -> float\n\nCross section for reaction MT at energy (MeV)."},
    {"xs_heating", (PyCFunction)AleaTHORNuclide_xs_heating, METH_VARARGS,
     "xs_heating(energy) -> float\n\nHeating number (MeV·barn) at energy."},
    {"heating_per_collision", (PyCFunction)AleaTHORNuclide_heating_per_collision, METH_VARARGS,
     "heating_per_collision(energy) -> float\n\nAverage energy deposited per collision (MeV)."},
    {"nu_bar", (PyCFunction)AleaTHORNuclide_nu_bar, METH_VARARGS,
     "nu_bar(energy) -> float\n\nAverage neutrons per fission at energy (MeV)."},
    {"doppler_broaden", (PyCFunction)AleaTHORNuclide_doppler_broaden, METH_VARARGS,
     "doppler_broaden(kT_target) -> None\n\n"
     "Doppler-broaden cross sections in-place to temperature kT_target (MeV).\n"
     "Only works on owned (non-cached) nuclides."},
    {"urr_factors", (PyCFunction)AleaTHORNuclide_urr_factors, METH_VARARGS,
     "urr_factors(energy, xi) -> dict or None\n\n"
     "Get URR probability table factors. Returns None if URR doesn't apply.\n"
     "Returns dict with keys: total, elastic, fission, capture, heating."},
    {"photon_xs_incoherent", (PyCFunction)AleaTHORNuclide_photon_xs_incoherent, METH_VARARGS,
     "photon_xs_incoherent(energy) -> float\n\nCompton (incoherent) scattering cross section (barns)."},
    {"photon_xs_coherent", (PyCFunction)AleaTHORNuclide_photon_xs_coherent, METH_VARARGS,
     "photon_xs_coherent(energy) -> float\n\nRayleigh (coherent) scattering cross section (barns)."},
    {"photon_xs_photoelectric", (PyCFunction)AleaTHORNuclide_photon_xs_photoelectric, METH_VARARGS,
     "photon_xs_photoelectric(energy) -> float\n\nPhotoelectric cross section (barns)."},
    {"photon_xs_pair", (PyCFunction)AleaTHORNuclide_photon_xs_pair, METH_VARARGS,
     "photon_xs_pair(energy) -> float\n\nPair production cross section (barns)."},
    {"energy_grid", (PyCFunction)AleaTHORNuclide_get_energy_grid, METH_NOARGS,
     "energy_grid() -> list of float\n\nGet the energy grid (MeV)."},
    {"reactions", (PyCFunction)AleaTHORNuclide_get_reactions, METH_NOARGS,
     "reactions() -> list of dict\n\nGet list of reactions with mt, q_value, ty."},
    {"reaction_yield", (PyCFunction)AleaTHORNuclide_reaction_yield, METH_VARARGS,
     "reaction_yield(mt, energy) -> float\n\nGet neutron yield for reaction MT at energy."},
    {NULL}
};

static PyTypeObject AleaTHORNuclideType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.Nuclide",
    .tp_doc = PyDoc_STR("Decoded ACE nuclide with cross-section lookup."),
    .tp_basicsize = sizeof(AleaTHORNuclideObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = AleaTHORNuclide_new,
    .tp_init = (initproc)AleaTHORNuclide_init,
    .tp_dealloc = (destructor)AleaTHORNuclide_dealloc,
    .tp_methods = AleaTHORNuclide_methods,
    .tp_getset = AleaTHORNuclide_getsetters,
};

/* ============================================================================
 * NucMaterial Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_material_t* mat;
    PyObject* refs;   /* List to keep nuclide refs alive */
} AleaTHORNucMaterialObject;

static void AleaTHORNucMaterial_dealloc(AleaTHORNucMaterialObject* self) {
    if (self->mat) {
        alea_nuc_material_destroy(self->mat);
    }
    Py_XDECREF(self->refs);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORNucMaterial_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    AleaTHORNucMaterialObject* self = (AleaTHORNucMaterialObject*)type->tp_alloc(type, 0);
    if (self) {
        self->mat = NULL;
        self->refs = NULL;
    }
    return (PyObject*)self;
}

static int AleaTHORNucMaterial_init(AleaTHORNucMaterialObject* self, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    self->mat = alea_nuc_material_create();
    if (!self->mat) {
        PyErr_SetString(PyExc_MemoryError, "Failed to create nuclear material");
        return -1;
    }
    self->refs = PyList_New(0);
    if (!self->refs) return -1;
    return 0;
}

static PyObject* AleaTHORNucMaterial_add(AleaTHORNucMaterialObject* self, PyObject* args) {
    AleaTHORNuclideObject* nuc_obj;
    double number_density;
    if (!PyArg_ParseTuple(args, "O!d", &AleaTHORNuclideType, &nuc_obj, &number_density))
        return NULL;

    if (!self->mat) {
        PyErr_SetString(PyExc_RuntimeError, "Material not initialized");
        return NULL;
    }
    if (!nuc_obj->nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized");
        return NULL;
    }

    alea_error_t err = alea_nuc_material_add(self->mat, nuc_obj->nuc, number_density);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_RuntimeError, "Failed to add nuclide: %s", alea_error_string(err));
        return NULL;
    }

    /* Keep nuclide alive */
    PyList_Append(self->refs, (PyObject*)nuc_obj);
    Py_RETURN_NONE;
}

/* Macroscopic cross sections */
static PyObject* AleaTHORNucMaterial_xs_total(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_total(self->mat, energy));
}

static PyObject* AleaTHORNucMaterial_xs_absorption(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_absorption(self->mat, energy));
}

static PyObject* AleaTHORNucMaterial_xs_elastic(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mat_xs_elastic(self->mat, energy));
}

static PyObject* AleaTHORNucMaterial_mean_free_path(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy;
    if (!PyArg_ParseTuple(args, "d", &energy)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mean_free_path(self->mat, energy));
}

static PyObject* AleaTHORNucMaterial_sample_distance(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_sample_distance(self->mat, energy, xi));
}

static PyObject* AleaTHORNucMaterial_sample_nuclide(AleaTHORNucMaterialObject* self, PyObject* args) {
    double energy, xi;
    if (!PyArg_ParseTuple(args, "dd", &energy, &xi)) return NULL;
    if (!self->mat) { PyErr_SetString(PyExc_RuntimeError, "Material not initialized"); return NULL; }

    alea_nuc_nuclide_t* out_nuc = NULL;
    int idx = alea_nuc_sample_nuclide(self->mat, energy, xi, &out_nuc);
    if (idx < 0 || !out_nuc) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to sample nuclide");
        return NULL;
    }

    return Py_BuildValue("(is)", idx, out_nuc->zaid);
}

static PyMethodDef AleaTHORNucMaterial_methods[] = {
    {"add", (PyCFunction)AleaTHORNucMaterial_add, METH_VARARGS,
     "add(nuclide, number_density)\n\nAdd a nuclide component with number density (atoms/barn-cm)."},
    {"xs_total", (PyCFunction)AleaTHORNucMaterial_xs_total, METH_VARARGS,
     "xs_total(energy) -> float\n\nMacroscopic total cross section (cm^-1)."},
    {"xs_absorption", (PyCFunction)AleaTHORNucMaterial_xs_absorption, METH_VARARGS,
     "xs_absorption(energy) -> float\n\nMacroscopic absorption cross section (cm^-1)."},
    {"xs_elastic", (PyCFunction)AleaTHORNucMaterial_xs_elastic, METH_VARARGS,
     "xs_elastic(energy) -> float\n\nMacroscopic elastic scattering cross section (cm^-1)."},
    {"mean_free_path", (PyCFunction)AleaTHORNucMaterial_mean_free_path, METH_VARARGS,
     "mean_free_path(energy) -> float\n\nMean free path (cm) at energy (MeV)."},
    {"sample_distance", (PyCFunction)AleaTHORNucMaterial_sample_distance, METH_VARARGS,
     "sample_distance(energy, xi) -> float\n\nSample distance to next collision (cm). xi in [0,1)."},
    {"sample_nuclide", (PyCFunction)AleaTHORNucMaterial_sample_nuclide, METH_VARARGS,
     "sample_nuclide(energy, xi) -> (index, zaid)\n\nSample which nuclide is hit. xi in [0,1)."},
    {NULL}
};

static PyTypeObject AleaTHORNucMaterialType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.NucMaterial",
    .tp_doc = PyDoc_STR("Nuclear material composition for transport calculations."),
    .tp_basicsize = sizeof(AleaTHORNucMaterialObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = AleaTHORNucMaterial_new,
    .tp_init = (initproc)AleaTHORNucMaterial_init,
    .tp_dealloc = (destructor)AleaTHORNucMaterial_dealloc,
    .tp_methods = AleaTHORNucMaterial_methods,
};

/* ============================================================================
 * Multigroup Python Type
 * ============================================================================ */

typedef struct {
    PyObject_HEAD
    alea_nuc_multigroup_t* mg;
    PyObject* spectrum_callable;   /* Python callable for custom spectrum, or NULL */
} AleaTHORMultigroupObject;

/* Trampoline: called from C, dispatches to the Python callable stored in ctx */
static double multigroup_spectrum_trampoline(double E, void* ctx) {
    PyObject* callable = (PyObject*)ctx;
    PyObject* result = PyObject_CallFunction(callable, "d", E);
    if (!result) {
        /* Can't propagate Python exceptions through C; return 0 and print */
        PyErr_Print();
        return 0.0;
    }
    double val = PyFloat_AsDouble(result);
    Py_DECREF(result);
    if (PyErr_Occurred()) {
        PyErr_Print();
        return 0.0;
    }
    return val;
}

static void AleaTHORMultigroup_dealloc(AleaTHORMultigroupObject* self) {
    if (self->mg) {
        alea_nuc_mg_destroy(self->mg);
    }
    Py_XDECREF(self->spectrum_callable);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* AleaTHORMultigroup_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    AleaTHORMultigroupObject* self = (AleaTHORMultigroupObject*)type->tp_alloc(type, 0);
    if (self) {
        self->mg = NULL;
        self->spectrum_callable = NULL;
    }
    return (PyObject*)self;
}

static int AleaTHORMultigroup_init(AleaTHORMultigroupObject* self, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {"bounds", NULL};
    PyObject* bounds_obj;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &bounds_obj))
        return -1;

    /* Convert bounds list to double array */
    PyObject* bounds_seq = PySequence_Fast(bounds_obj, "bounds must be a sequence");
    if (!bounds_seq) return -1;

    Py_ssize_t n = PySequence_Fast_GET_SIZE(bounds_seq);
    if (n < 2) {
        Py_DECREF(bounds_seq);
        PyErr_SetString(PyExc_ValueError, "Need at least 2 group boundaries");
        return -1;
    }

    double* bounds = (double*)malloc(sizeof(double) * n);
    if (!bounds) {
        Py_DECREF(bounds_seq);
        PyErr_NoMemory();
        return -1;
    }

    for (Py_ssize_t i = 0; i < n; i++) {
        bounds[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bounds_seq, i));
        if (PyErr_Occurred()) {
            free(bounds);
            Py_DECREF(bounds_seq);
            return -1;
        }
    }
    Py_DECREF(bounds_seq);

    int n_groups = (int)(n - 1);
    self->mg = alea_nuc_mg_create(n_groups, bounds);
    free(bounds);

    if (!self->mg) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create multigroup structure");
        return -1;
    }

    return 0;
}

static PyObject* AleaTHORMultigroup_collapse(AleaTHORMultigroupObject* self, PyObject* args) {
    AleaTHORNuclideObject* nuc_obj;
    if (!PyArg_ParseTuple(args, "O!", &AleaTHORNuclideType, &nuc_obj)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    if (!nuc_obj->nuc) { PyErr_SetString(PyExc_RuntimeError, "Nuclide not initialized"); return NULL; }

    alea_error_t err = alea_nuc_mg_collapse(self->mg, nuc_obj->nuc);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_RuntimeError, "Multigroup collapse failed: %s", alea_error_string(err));
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject* AleaTHORMultigroup_scatter(AleaTHORMultigroupObject* self, PyObject* args) {
    int g_from, g_to;
    if (!PyArg_ParseTuple(args, "ii", &g_from, &g_to)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mg_scatter(self->mg, g_from, g_to));
}

static PyObject* AleaTHORMultigroup_scatter_adjoint(AleaTHORMultigroupObject* self, PyObject* args) {
    int g_from, g_to;
    if (!PyArg_ParseTuple(args, "ii", &g_from, &g_to)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyFloat_FromDouble(alea_nuc_mg_scatter_adjoint(self->mg, g_from, g_to));
}

static PyObject* AleaTHORMultigroup_sample_scatter(AleaTHORMultigroupObject* self, PyObject* args) {
    int g_from, adjoint = 0;
    double xi;
    if (!PyArg_ParseTuple(args, "id|p", &g_from, &xi, &adjoint)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyLong_FromLong(alea_nuc_mg_sample_scatter(self->mg, g_from, xi, adjoint));
}

static PyObject* AleaTHORMultigroup_set_spectrum(AleaTHORMultigroupObject* self, PyObject* args) {
    PyObject* callable;
    if (!PyArg_ParseTuple(args, "O", &callable)) return NULL;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }

    if (callable == Py_None) {
        /* Reset to default spectrum */
        alea_nuc_mg_set_spectrum(self->mg, NULL, NULL);
        Py_CLEAR(self->spectrum_callable);
    } else {
        if (!PyCallable_Check(callable)) {
            PyErr_SetString(PyExc_TypeError, "spectrum must be callable or None");
            return NULL;
        }
        Py_INCREF(callable);
        Py_XDECREF(self->spectrum_callable);
        self->spectrum_callable = callable;
        alea_nuc_mg_set_spectrum(self->mg, multigroup_spectrum_trampoline, callable);
    }
    Py_RETURN_NONE;
}

/* Get group constants as dict */
static PyObject* AleaTHORMultigroup_get_data(AleaTHORMultigroupObject* self, PyObject* args) {
    (void)args;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }

    int ng = self->mg->n_groups;
    PyObject* result = PyDict_New();
    if (!result) return NULL;

    /* Helper macro to create list from double array */
    #define MG_LIST(name, arr, size) do { \
        PyObject* list = PyList_New(size); \
        if (!list) { Py_DECREF(result); return NULL; } \
        for (int i = 0; i < (size); i++) \
            PyList_SET_ITEM(list, i, PyFloat_FromDouble((arr)[i])); \
        PyDict_SetItemString(result, name, list); \
        Py_DECREF(list); \
    } while(0)

    MG_LIST("bounds", self->mg->bounds, ng + 1);
    MG_LIST("sigma_t", self->mg->sigma_t, ng);
    MG_LIST("sigma_a", self->mg->sigma_a, ng);
    MG_LIST("sigma_s", self->mg->sigma_s, ng);
    MG_LIST("sigma_f", self->mg->sigma_f, ng);
    MG_LIST("nu_sigma_f", self->mg->nu_sigma_f, ng);
    MG_LIST("chi", self->mg->chi, ng);

    /* Scattering matrix as list of lists */
    PyObject* scatter_mat = PyList_New(ng);
    if (!scatter_mat) { Py_DECREF(result); return NULL; }
    for (int g = 0; g < ng; g++) {
        PyObject* row = PyList_New(ng);
        if (!row) { Py_DECREF(scatter_mat); Py_DECREF(result); return NULL; }
        for (int gp = 0; gp < ng; gp++) {
            PyList_SET_ITEM(row, gp, PyFloat_FromDouble(self->mg->scatter[g * ng + gp]));
        }
        PyList_SET_ITEM(scatter_mat, g, row);
    }
    PyDict_SetItemString(result, "scatter_matrix", scatter_mat);
    Py_DECREF(scatter_mat);

    #undef MG_LIST

    PyDict_SetItemString(result, "n_groups", PyLong_FromLong(ng));

    return result;
}

static PyObject* AleaTHORMultigroup_get_n_groups(AleaTHORMultigroupObject* self, void* closure) {
    (void)closure;
    if (!self->mg) { PyErr_SetString(PyExc_RuntimeError, "Multigroup not initialized"); return NULL; }
    return PyLong_FromLong(self->mg->n_groups);
}

static PyGetSetDef AleaTHORMultigroup_getsetters[] = {
    {"n_groups", (getter)AleaTHORMultigroup_get_n_groups, NULL, "Number of energy groups", NULL},
    {NULL}
};

static PyMethodDef AleaTHORMultigroup_methods[] = {
    {"collapse", (PyCFunction)AleaTHORMultigroup_collapse, METH_VARARGS,
     "collapse(nuclide)\n\nCollapse continuous-energy cross sections into multigroup constants."},
    {"set_spectrum", (PyCFunction)AleaTHORMultigroup_set_spectrum, METH_VARARGS,
     "set_spectrum(fn)\n\n"
     "Set a custom weighting spectrum for group collapse.\n"
     "fn must be a callable taking energy (MeV) and returning the spectrum value,\n"
     "or None to reset to the default (Maxwellian + 1/E + fission)."},
    {"scatter", (PyCFunction)AleaTHORMultigroup_scatter, METH_VARARGS,
     "scatter(g_from, g_to) -> float\n\nForward scattering matrix element."},
    {"scatter_adjoint", (PyCFunction)AleaTHORMultigroup_scatter_adjoint, METH_VARARGS,
     "scatter_adjoint(g_from, g_to) -> float\n\nAdjoint scattering matrix element."},
    {"sample_scatter", (PyCFunction)AleaTHORMultigroup_sample_scatter, METH_VARARGS,
     "sample_scatter(g_from, xi, adjoint=False) -> int\n\nSample outgoing group from scattering."},
    {"get_data", (PyCFunction)AleaTHORMultigroup_get_data, METH_NOARGS,
     "get_data() -> dict\n\nGet all multigroup constants as a dictionary."},
    {NULL}
};

static PyTypeObject AleaTHORMultigroupType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_alea.Multigroup",
    .tp_doc = PyDoc_STR("Multigroup cross sections and scattering matrix."),
    .tp_basicsize = sizeof(AleaTHORMultigroupObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = AleaTHORMultigroup_new,
    .tp_init = (initproc)AleaTHORMultigroup_init,
    .tp_dealloc = (destructor)AleaTHORMultigroup_dealloc,
    .tp_methods = AleaTHORMultigroup_methods,
    .tp_getset = AleaTHORMultigroup_getsetters,
};

/* ============================================================================
 * Module-level nucdata functions
 * ============================================================================ */

static PyObject* mod_parse_zaid(PyObject* self, PyObject* args) {
    (void)self;
    const char* zaid;
    if (!PyArg_ParseTuple(args, "s", &zaid)) return NULL;

    int Z, A, meta;
    alea_nuc_table_type_t type;
    alea_error_t err = alea_nuc_parse_zaid(zaid, &Z, &A, &meta, &type);
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_ValueError, "Invalid ZAID '%s': %s", zaid, alea_error_string(err));
        return NULL;
    }

    const char* type_str;
    switch (type) {
        case ALEA_NUC_TABLE_CONTINUOUS_NEUTRON: type_str = "continuous_neutron"; break;
        case ALEA_NUC_TABLE_PHOTOATOMIC:        type_str = "photoatomic"; break;
        case ALEA_NUC_TABLE_PHOTONUCLEAR:       type_str = "photonuclear"; break;
        case ALEA_NUC_TABLE_THERMAL_SAB:        type_str = "thermal_sab"; break;
        case ALEA_NUC_TABLE_ELECTRON:           type_str = "electron"; break;
        default:                                 type_str = "unknown"; break;
    }

    return Py_BuildValue("{s:i, s:i, s:i, s:s}",
        "Z", Z, "A", A, "metastable", meta, "type", type_str);
}

static PyObject* mod_reaction_classify(PyObject* self, PyObject* args) {
    (void)self;
    int mt;
    if (!PyArg_ParseTuple(args, "i", &mt)) return NULL;

    alea_nuc_reaction_class_t cls = alea_nuc_reaction_classify(mt);
    const char* cls_str;
    switch (cls) {
        case ALEA_NUC_RXN_ABSORPTION: cls_str = "absorption"; break;
        case ALEA_NUC_RXN_SCATTER:    cls_str = "scatter"; break;
        case ALEA_NUC_RXN_MULTIPLY:   cls_str = "multiply"; break;
        default:                       cls_str = "unknown"; break;
    }

    return PyUnicode_FromString(cls_str);
}
