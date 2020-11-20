#ifndef TOPOPTLIB_HPP
#define TOPOPTLIB_HPP

#include <Python.h>
#include <vector>
#include <string>
#include <petsc.h>

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
        // data storage standard var
        double xc_w[11] = {0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        int nxyz_w[3] = {65, 33, 33};
        double Emin_w  = 0.0;
        double Emax_w = 1.0;
        double nu_w = 0.3;
        double penal_w = 3.0;
        int maxIter_w = 40;
        int filter_w = 1;
        double rmin_w = 0.08;
        double volumefrac_w = 0.12;

        // not needed
        double b_w[6];

        // continuation, projections
        int continuation_w = 0;
        int projection_w = 0;
        int localVolume_w = 0;
        double betaInit_w = 1.0;
        double betaFinal_w = 64.0;
        double eta_w = 0.5;
        double penalfinal_w;
        double stepsize_w;

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
        int test_w = 0;

        // xPassive
        std::vector<double> xPassive_w;
        int nael; // number of active design variables
        int nrel; // number of rigid elements
        int nsel; // number of solid elements

        // constraints
        int m = 1;
        double Rlocvol_w, alpha_w;

        // Loadcases
        int nL = 1;
        std::vector<std::vector<BC>> loadcases_list;
        
        // parametrization
        PyObject *para_func = NULL;    

        // test function
        PyObject *test_func = NULL;

        // Callback function
        PyObject *obj_func = NULL;
        PyObject *obj_sens_func = NULL;
        PyObject *const_func = NULL;
        PyObject *const_sens_func = NULL;

        // 
        void updatecounts() {
            
            int acount = 0;
            int scount = 0;
            int rcount = 0;
            // count the number of active, solid and rigid elements
            // active is -1, passive is 1, solid is 2, rigid is 3
            for (unsigned i = 0; i < xPassive_w.size(); i++) {
                if (xPassive_w[i] == -1.0) {
                    acount++;
                }
                if (xPassive_w[i] == 2.0) {
                    scount++;
                } 
                if (xPassive_w[i] == 3.0) {
                    rcount++;
                } 
            }
            
            // reset
            nael = acount;
            nsel = scount;
            nrel = rcount;
            
        }

        void check_ev(double trueFx, double runtime, double memory) {
            
            //printf("el: %i\n", el);
            //printf("uKu: %f\n", uKu);

            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("ddd", trueFx, runtime, memory);
            
            // PyEval_CallObject
            PyEval_CallObject(test_func, arglist);

            // Parse
            //double result = PyFloat_AsDouble(results);

            //Py_XDECREF(results);
            Py_DECREF(arglist);

            //return;
        }

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

        double obj_ev(double comp, double sumXP) {
            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("dd", comp, sumXP);
            
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

        double const_ev(double comp, double sumXP) {
            // Py_BuildValue
            PyObject *arglist = Py_BuildValue("dd", comp, sumXP);
            
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

int solve(DataObj data);
#endif