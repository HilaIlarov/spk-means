
#define PY_SSIZE_T_CLEAN /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>      /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                        <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <math.h>        /* include <Python.h> has to be before any standard headers are included */
#include "spkmeans.h"
#include <stdio.h>

static PyObject *c_to_py(double **data_points, int N, int dim)
{
    PyObject *python_data_points;
    PyObject *python_vector;
    int i;
    int j;
    python_data_points = PyList_New(N);
    for (i = 0; i < N; i++)
    {
        python_vector = PyList_New(dim);
        for (j = 0; j < dim; j++)
        {
            PyList_SetItem(python_vector, j, PyFloat_FromDouble(data_points[i][j]));
        }
        PyList_SetItem(python_data_points, i, python_vector);
    }
    return python_data_points;
}

double **py_to_c(PyObject *data_points, int N, int dim)
{
    double **c_data_points;
    PyObject *x_i;
    PyObject *x_i_j;
    int i;
    int j;
    c_data_points = callocation(N, dim + 1);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < dim; j++)
        {
            x_i = PyList_GetItem(data_points, i);
            x_i_j = PyList_GetItem(x_i, j);
            c_data_points[i][j] = PyFloat_AsDouble(x_i_j);
        }
    }
    return c_data_points;
}

static PyObject *send_T_to_python(PyObject *self, PyObject *args)
{
    int k;
    char *goal;
    int dim;
    int N;
    PyObject *vectors_py;
    PyObject *t_datapoints_py;
    double **vectors_c;
    double **t_datapoints_c;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if (!PyArg_ParseTuple(args, "isiiO", &k, &goal, &N, &dim, &vectors_py))
    {
        return NULL; /* In the CPython API, a NULL value is never valid for a
        PyObject* so it is used to signal that an error has occurred. */
    }
    vectors_c = py_to_c(vectors_py, N, dim);
    /* This builds the answer ("d" = Convert a C double to a Python floating point number) back into a python object */
    t_datapoints_c = spk(k, goal, vectors_c, N, dim, 0);
    k = (int)t_datapoints_c[N][0];
    t_datapoints_py = c_to_py(t_datapoints_c, N, k);
    free_array(vectors_c, N);
    free_array(t_datapoints_c, N);
    return Py_BuildValue("O", t_datapoints_py); /*  Py_BuildValue(...) returns a PyObject*  */
}

static PyObject *fit(PyObject *self, PyObject *args)
{
    /*Initializing*/
    int K;
    int max_iter;
    double *v_zero = NULL;
    double **x = NULL;
    double **mu = NULL;
    cluster *s = NULL;
    PyObject *x_py;
    PyObject *mu_py;
    PyObject *python_mu_list;
    int N = 0;
    int d = 0;
    int i = 1;

    if (!PyArg_ParseTuple(args, "iiiiOO", &N, &d, &K, &max_iter, &x_py, &mu_py))
    {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    /*x_array allocation size of N x d */
    /*x_array filling with vectors from the pyton list*/
    /*x[i][d] = j <=> x is in Sj  ;  init with -1*/
    x = py_to_c(x_py, N, d);
    for (i = 0; i < N; i++)
    {
        x[i][d] = -1;
    }

    /*allocating mus_array*/
    mu = py_to_c(mu_py, K, d);
    assert((K < N) && "k must be lower than N");
    v_zero = build_v_zero(v_zero, d);

    /*allocating S*/
    s = allocate_s(s, v_zero, K);

    /*calc mu*/
    mu = build_mu(x, mu, s, v_zero, N, K, d);
    /*transform from double ** to pylist*/
    python_mu_list = c_to_py(mu, K, d);

    /*free them all*/
    /*free x, mu, s, v_zero*/
    free_array(x, N);
    free_array(mu, K);
    free_s(s, K);
    free(v_zero);

    return Py_BuildValue("O", python_mu_list);
}

/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef Methods[] = {
    {"send_T_to_python",                                     /* the Python method name that will be used */
     (PyCFunction)send_T_to_python,                          /* the C-function that implements the Python function and returns static PyObject*  */
     METH_VARARGS,                                           /* flags indicating parameters accepted for this function */
     PyDoc_STR("run spk and send T matrix back to python")}, /*  The docstring for the function */
    {"fit",                                                  /* the Python method name that will be used */
     (PyCFunction)fit,                                       /* the C-function that implements the Python function and returns static PyObject*  */
     METH_VARARGS,                                           /* flags indicating parameters accepted for this function */
     PyDoc_STR("run hw2")},                                  /* The docstring for the function */
    {NULL, NULL, 0, NULL}                                    /* The last entry must be all NULL as shown to act as a sentinel. 
                                                                Python looks for this entry to know that all
                                                                of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans", /* name of module */
    NULL,       /* module documentation, may be NULL */
    -1,         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    Methods     /* the PyMethodDef array from before containing the methods of the extension */
};

/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
    {
        return NULL;
    }
    return m;
}