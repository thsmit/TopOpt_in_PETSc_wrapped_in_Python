#ifndef TOPOPTLIB_HPP
#define TOPOPTLIB_HPP

#include <Python.h>
#include <vector>

/*
Author: Thijs Smit, May 2020
Copyright (C) 2020 ETH Zurich

Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

class BC {
    public:

        int BCtype;
        int Para;

        std::vector<int> Checker_vec;
        std::vector<int> Setter_dof_vec;
        std::vector<double> Setter_val_vec;
};

// Creating class with C variables with default values
// check if class works
struct DataObj {
    // check PyObject_class?
    PyObject_HEAD
    
    public:
        // data storage
        double xc_w[9];
        double b_w[6];
        int nxyz_w[3] = {65, 33, 33};
        double Emin_w = 1.0e-9;
        double Emax_w = 1.0;
        double nu_w = 0.3;
        double penal_w = 3.0;
        int maxIter_w = 400;
        int filter_w = 1;
        double rmin_w = 0.08;
        double volumefrac_w = 0.12;

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

        // xPassive
        std::vector<double> xPassive_w;
        int nael; // number of active design variables

        // Loadcases
        int nL = 1;
        std::vector<std::vector<BC>> loadcases_list;
        
        // parametrization
        PyObject *para_func = NULL;    

        // parametrization
        PyObject *passive_func = NULL;

        // Callback function
        PyObject *obj_func = NULL;
        PyObject *obj_sens_func = NULL;
        PyObject *const_func = NULL;
        PyObject *const_sens_func = NULL;

        double para_ev(double lcx, double lcy, double lcz) {
            
            //printf("el: %i\n", el);
            //printf("uKu: %f\n", uKu);

            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("ddd", lcx, lcy, lcz);
            
            // PyEval_CallObject
            PyObject *results = PyEval_CallObject(para_func, arglist);

            // Parse
            double result = PyFloat_AsDouble(results);

            Py_XDECREF(results);
            Py_DECREF(arglist);

            return result;
        }

        double passive_ev(double el) {
            
            //printf("el: %i\n", el);
            //printf("uKu: %f\n", uKu);

            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("d", el);
            
            // PyEval_CallObject
            PyObject *results = PyEval_CallObject(passive_func, arglist);
            
            // Parse
            double result = PyFloat_AsDouble(results);

            Py_XDECREF(results);
            Py_DECREF(arglist);

            return result;
        }

        double obj_ev(double comp, double sumXP, double xp, double uKu) {
            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("dddd", comp, sumXP, xp, uKu);
            
            // PyEval_CallObject
            PyObject *results = PyEval_CallObject(obj_func, arglist);
            
            // Parse
            double result = PyFloat_AsDouble(results);

            Py_XDECREF(results);
            Py_DECREF(arglist);

            return result;
        }

        double obj_sens_ev(double sumXP, double xp, double uKu) {
            
            //printf("xp: %f\n", xp);
            //printf("uKu: %f\n", uKu);

            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("ddd", sumXP, xp, uKu);
            
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

        double const_ev(double comp, double sumXP, double xp, double uKu) {
            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("dddd", comp, sumXP, xp, uKu);
            
            // PyEval_CallObject
            PyObject *results = PyEval_CallObject(const_func, arglist);
            
            // Parse
            double result = PyFloat_AsDouble(results);

            Py_XDECREF(results);
            Py_DECREF(arglist);

            return result;
        }

        double const_sens_ev(double sumXP, double xp, double uKu) {
            
            //printf("xp: %f\n", xp);
            //printf("uKu: %f\n", uKu);

            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("ddd", sumXP, xp, uKu);
            
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

int solve(DataObj data);
#endif