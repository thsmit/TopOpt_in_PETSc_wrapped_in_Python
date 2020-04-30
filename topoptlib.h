#ifndef TOPOPTLIB_HPP
#define TOPOPTLIB_HPP

#include <Python.h>

/*
Author: Thijs Smit, April 2020
Copyright (C) 2020 ETH Zurich

Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

// Creating class with C variables with default values
struct DataObj {
    PyObject_HEAD
    
    // data storage
    double xc_w[6] = {0.0, 2.0, 0.0, 1.0, 0.0, 1.0};
    int nxyz_w[3] = {65, 33, 33};
    double Emin_w = 1.0e-9;
    double Emax_w = 1.0;
    double nu_w = 0.3;
    double penal_w = 3.0;
    int maxIter_w = 400;
    int filter_w = 1;
    double rmin_w = 0.08;

    // needed as members to be used in python script
    int nNodes;
    int nElements;
    int nDOF;

    int it;
    int nCores;
    int avMem;
    int peakMem;
    int cpuTime;
    double trueFx;
    double scaledFx;

    // Lock
    //PyThread_type_lock lock;

    // BC
    PyObject *loadcases_list = NULL;
    int nL = 1;

    // Passive elements
    //bool passive = false;
    PyObject *passive_func = NULL;

    // Callback function
    PyObject *obj_func = NULL;
    PyObject *obj_sens_func = NULL;
    PyObject *const_func = NULL;
    PyObject *const_sens_func = NULL;

    double passive_ev(int el) {
        
        //printf("el: %i\n", el);
        //printf("uKu: %f\n", uKu);

        // Py_BuildValue
        PyObject *arglist = Py_BuildValue("i", el);
        
        // PyEval_CallObject
        PyObject *results = PyEval_CallObject(passive_func, arglist);
        
        // Parse
        //double result;

        //result = PyFloat_AsDouble(results);
        
        double result = PyFloat_AsDouble(results);

        //printf("result: %f\n", result);
        Py_XDECREF(results);
        Py_DECREF(arglist);

        return result;
    }

    double obj_ev(double xp, double uKu) {
        // Py_BuildValue
        PyObject *arglist = Py_BuildValue("dd", xp, uKu);
        
        // PyEval_CallObject
        PyObject *results = PyEval_CallObject(obj_func, arglist);
        
        // Parse
        double result = PyFloat_AsDouble(results);

        Py_XDECREF(results);
        Py_DECREF(arglist);

        return result;
    }

    double obj_sens_ev(double xp, double uKu) {
        
        //printf("xp: %f\n", xp);
        //printf("uKu: %f\n", uKu);

        // Py_BuildValue
        PyObject *arglist = Py_BuildValue("dd", xp, uKu);
        
        // PyEval_CallObject
        PyObject *results = PyEval_CallObject(obj_sens_func, arglist);
        
        // Parse
        //double result;

        //result = PyFloat_AsDouble(results);
        
        double result = PyFloat_AsDouble(results);

        //printf("result: %f\n", result);
        Py_XDECREF(results);
        Py_DECREF(arglist);

        return result;
    }

     double const_ev(double xp, double uKu) {
        // Py_BuildValue
        PyObject *arglist = Py_BuildValue("dd", xp, uKu);
        
        // PyEval_CallObject
        PyObject *results = PyEval_CallObject(const_func, arglist);
        
        // Parse
        double result = PyFloat_AsDouble(results);

        Py_XDECREF(results);
        Py_DECREF(arglist);

        return result;
    }

    double const_sens_ev(double xp, double uKu) {
        
        //printf("xp: %f\n", xp);
        //printf("uKu: %f\n", uKu);

        // Py_BuildValue
        PyObject *arglist = Py_BuildValue("dd", xp, uKu);
        
        // PyEval_CallObject
        PyObject *results = PyEval_CallObject(const_sens_func, arglist);
        
        // Parse
        //double result;

        //result = PyFloat_AsDouble(results);
        
        double result = PyFloat_AsDouble(results);

        //printf("result: %f\n", result);
        Py_XDECREF(results);
        Py_DECREF(arglist);

        return result;
    }

};



/*
class Data {
    public:

        // default input, is overwritten by user input data
        double xc_w[6] = {0.0, 2.0, 0.0, 1.0, 0.0, 1.0};
        int nxyz_w[3] = {65, 33, 33};
        double Emin_w = 1.0e-9;
        double Emax_w = 1.0;
        double nu_w = 0.3;
        double penal_w = 3.0;
        int nLC = 1;
        int maxIter_w = 400;

        // Callback functions
        //int obj(xp, uKu) {
        //    return NULL;
        //}    
};
*/

int solve(DataObj data);
#endif