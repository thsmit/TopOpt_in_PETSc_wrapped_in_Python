#include "topoptlib.h" // for wrapper
#include "Filter.h"
#include "LocalVolume.h"
#include "LinearElasticity.h"
#include "MMA.h"
#include "MPIIO.h"
#include "TopOpt.h"
#include "mpi.h"
#include <petsc.h>

/*
Modified: Thijs Smit, Dec 2020

Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

Updated: June 2019, Niels Aage
Copyright (C) 2013-2019,

Disclaimer:
The authors reserves all rights but does not guaranty that the code is
free from errors. Furthermore, we shall not be liable in any event
caused by the use of the program.
*/

static char help[] = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";

int solve(DataObj data) {

    // HACK to fake zero run time input data...
    int argc = 0;
    char **argv = 0;

    // Error code for debugging
    PetscErrorCode ierr;

    // Initialize PETSc / MPI and pass input arguments to PETSc
    PetscInitialize(&argc, &argv, PETSC_NULL, help);

    // Capture runtime
    double rt1, rt2;

    // Monitor memory usage
    PetscMemorySetGetMaximumUsage();

    // start timer total runtime
    rt1 = MPI_Wtime();

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
    TopOpt* opt = new TopOpt(data);

    // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes, data);

    // STEP 3: THE FILTERING
    //Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->xPassive, opt->filter, opt->rmin);

    // Initialize local volume constraint
    LocalVolume* local;
    if (opt->localVolumeStatus) {
        local = new LocalVolume(opt->da_nodes, opt->x, opt->xPassive, data);
    }

    MPIIO* output;
    if (opt->robustStatus) {
        // STEP 4: VISUALIZATION USING VTK
        output = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 5, "x, xTilde, xPhysEro, xPhys, xPhysDil");
    } else {
        output = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 3, "x, xTilde, xPhys");
    }

    // STEP 5: THE OPTIMIZER MMA
    MMA* mma;
    PetscInt itr = 0;
    opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !

    if (opt->robustStatus) {
        // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
        ierr = filter->FilterProjectRobust(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta);
        CHKERRQ(ierr);
    } else {
        ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta);
        CHKERRQ(ierr);
    }

    // Set variables in xPhys, xPhysEro, xPhysDil
    opt->SetVariables(opt->xPhys, opt->xPassive);
    if (opt->robustStatus) {
        opt->SetVariables(opt->xPhysEro, opt->xPassive);
        opt->SetVariables(opt->xPhysDil, opt->xPassive);
    }

    if (opt->robustStatus) {
        // print initial condition to vtk
        output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, itr);
    } else {
        // print initial condition to vtk
        output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr);
    }


    // STEP 7: OPTIMIZATION LOOP
    PetscScalar ch = 1.0;
    double      t1, t2;
    //while (ch > opt->tol || opt->penal < opt->penalFin) {
    //while ((itr < opt->maxItr && ch > opt->tol) || opt->penal < opt->penalFin) {
    while (itr < opt->maxItr && ch > opt->tol) {
    //while (itr < opt->maxItr && ch > opt->tol || opt->beta < opt->betaFinal) {
        // Update iteration counter
        itr++;

		if (opt->continuationStatus && itr % opt->IterProj == 0)
        //if (opt->continuationStatus && ((opt->penal != 1.0 && itr % opt->IterProj == 0) || (opt->penal == 1.0 && itr > 150)))
        //if (opt->continuationStatus && ch < opt->tol)
        {
			opt->penal = PetscMin(opt->penal + opt->penalStep, opt->penalFin);
			PetscPrintf(PETSC_COMM_WORLD,"It.: %i, Beta: %3.2f, Penal: %2.2f   movlim: %f\n", itr,opt->beta,opt->penal,opt->movlim);
		}

        // start timer
        t1 = MPI_Wtime();

        if (opt->robustStatus) {
            // Scaling the upper bound volume constraint of the Dilated Design
            if (itr % 20 == 0){
                PetscScalar VolDil = 0.0;
                PetscScalar VolBlue = 0.0;
                VecSum(opt->xPhysDil, &VolDil);
                VecSum(opt->xPhys, &VolBlue);
                opt->volfrac = ((VolDil - data.nrel * 10.0)  / (VolBlue - data.nrel * 10.0) ) * opt->volfracREF;
                PetscPrintf(PETSC_COMM_WORLD,"volfrac: %f\n", opt->volfrac);
            }
            // Compute obj+const+sens
            ierr = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), &(opt->gx[0]), opt->dfdx, opt->dgdx[0],
                                                                    opt->xPhysEro, opt->xPhysDil, opt->Emin, opt->Emax, opt->penal,
                                                                    opt->volfrac, data);
            CHKERRQ(ierr);
        } else {
            // Compute obj+const+sens
            ierr = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), &(opt->gx[0]), opt->dfdx, opt->dgdx[0],
                                                                 opt->xPhys, opt->xPhys, opt->Emin, opt->Emax, opt->penal,
                                                                 opt->volfrac, data);
            CHKERRQ(ierr);
        }

        // Compute objective scale
        if (itr == 1 || (itr % opt->IterProj == 0 && opt->continuationStatus)) {
            opt->fscale = 10.0 / opt->fx;
        }
        // Scale objectie and sens
        opt->fx = opt->fx * opt->fscale;
        VecScale(opt->dfdx, opt->fscale);

        if (opt->robustStatus) {
            // Filter sensitivities (chainrule)
            ierr = filter->GradientsRobust(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta,
                                    opt->eta);
            CHKERRQ(ierr);
        } else {
            ierr = filter->Gradients(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta,
                                    opt->eta);
            CHKERRQ(ierr);
        }

        // Calculate g and dgdx for the local volume constraint
        if (opt->localVolumeStatus) {
            ierr = local->Constraint(opt->x, &(opt->gx[1]), opt->dgdx[1]);
            CHKERRQ(ierr);
        }

        if (opt->xPassiveStatus) {
            // map vectors
            opt->UpdateVariables(1, opt->x, opt->xMMA);
            opt->UpdateVariables(1, opt->dfdx, opt->dfdxMMA);
            for (PetscInt i = 0; i < opt->m; i++) {
                opt->UpdateVariables(1, opt->dgdx[i], opt->dgdxMMA[i]);
            }
            opt->UpdateVariables(1, opt->xmin, opt->xminMMA);
            opt->UpdateVariables(1, opt->xmax, opt->xmaxMMA);
            opt->UpdateVariables(1, opt->xold, opt->xoldMMA);

            // Sets outer movelimits on design variables
            ierr = mma->SetOuterMovelimit(opt->Xmin, opt->Xmax, opt->movlim, opt->xMMA, opt->xminMMA, opt->xmaxMMA);
            CHKERRQ(ierr);

            // Update design by MMA
            ierr = mma->Update(opt->xMMA, opt->dfdxMMA, opt->gx, opt->dgdxMMA, opt->xminMMA, opt->xmaxMMA);
            CHKERRQ(ierr);

            // Inf norm on the design change
            ch = mma->DesignChange(opt->xMMA, opt->xoldMMA);

            opt->UpdateVariables(-1, opt->x, opt->xMMA);
            opt->UpdateVariables(-1, opt->dfdx, opt->dfdxMMA);
            for (PetscInt i = 0; i < opt->m; i++) {
                opt->UpdateVariables(-1, opt->dgdx[i], opt->dgdxMMA[i]);
            }
            opt->UpdateVariables(-1, opt->xmin, opt->xminMMA);
            opt->UpdateVariables(-1, opt->xmax, opt->xmaxMMA);
            opt->UpdateVariables(-1, opt->xold, opt->xoldMMA);

        } else {

            // Sets outer movelimits on design variables
            ierr = mma->SetOuterMovelimit(opt->Xmin, opt->Xmax, opt->movlim, opt->x, opt->xmin, opt->xmax);
            CHKERRQ(ierr);

            // Update design by MMA
            ierr = mma->Update(opt->x, opt->dfdx, opt->gx, opt->dgdx, opt->xmin, opt->xmax);
            CHKERRQ(ierr);

            // Inf norm on the design change
            ch = mma->DesignChange(opt->x, opt->xold);

        }

        // Increase beta if needed
        PetscBool changeBeta = PETSC_FALSE;
        if (opt->projectionFilter) {
        //if (opt->projectionFilter && opt->penal == opt->penalFin) {
        //if (opt->projectionFilter && itr > 150) {
            changeBeta = filter->IncreaseBeta(&(opt->beta), opt->betaFinal, opt->gx[0], itr, ch);
        }

        if (opt->robustStatus) {
            // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
            ierr = filter->FilterProjectRobust(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta);
            CHKERRQ(ierr);
        } else {
            ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, opt->projectionFilter, opt->beta, opt->eta);
            CHKERRQ(ierr);
        }

        // Set variables in xPhys, xPhysEro, xPhysDil
        opt->SetVariables(opt->xPhys, opt->xPassive);
        if (opt->robustStatus) {
            opt->SetVariables(opt->xPhysEro, opt->xPassive);
            opt->SetVariables(opt->xPhysDil, opt->xPassive);
        }
        // Discreteness measure
        PetscScalar mnd = filter->GetMND(opt->xPhys);

        // stop timer
        t2 = MPI_Wtime();

        // Print to screen
        PetscPrintf(PETSC_COMM_WORLD,
                    "It.: %i, True fx: %f, Scaled fx: %f, gx[0]: %f, ch.: %f, "
                    "mnd.: %f, time: %f\n",
                    itr, opt->fx / opt->fscale, opt->fx, opt->gx[0], ch, mnd, t2 - t1);

        // Write field data: first 10 iterations and then every 20th
        if (itr < 11 || itr % 20 == 0 || changeBeta) {
            if (opt->robustStatus) {
                // print initial condition to vtk
                output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, itr);
            } else {
                // print initial condition to vtk
                output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr);
            }
        }



        // Dump data needed for restarting code at termination
        //if (itr % 10 == 0) {
            //opt->WriteRestartFiles(&itr, mma);
            //physics->WriteRestartFiles();
        //}
    }

    // Write restart WriteRestartFiles
    //opt->WriteRestartFiles(&itr, mma);
    //physics->WriteRestartFiles();

    // Dump final design
    output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhysEro, opt->xPhys, opt->xPhysDil, itr + 1);

    // stop timer total runtime
    rt2 = MPI_Wtime();

    // Extract memory use
    PetscLogDouble mem;
    PetscMemoryGetMaximumUsage(&mem);

    // Call test function to send data to wrapper
    if (opt->testStatus) {
        data.check_ev(opt->fx / opt->fscale, rt2 - rt1, mem);
    }

    //// Print output
    PetscScalar xPhys_sum;
    VecSum(opt->xPhys, &(xPhys_sum));
    //PetscPrintf(PETSC_COMM_WORLD, "# xPhys_sum: %f\n", xPhys_sum);
    xPhys_sum = xPhys_sum - data.nrel * 10.0 - data.nsel * 1.0;
    //PetscPrintf(PETSC_COMM_WORLD, "# xPhys_sum1: %f\n", xPhys_sum);

    // Calculate porosity
    //PetscScalar Poro = 1.0 - (xPhys_sum / ((opt->nxyz[0] - 1) * (opt->nxyz[1] - 1) * (opt->nxyz[2] - 1)));
    PetscScalar Poro = 1.0 - (xPhys_sum / data.nael);

    PetscPrintf(PETSC_COMM_WORLD, "######################## Final output ########################\n");
    PetscPrintf(PETSC_COMM_WORLD, "# Porosity: %f\n", Poro);
    PetscPrintf(PETSC_COMM_WORLD, "# Total Volume: %f\n", xPhys_sum);
    PetscPrintf(PETSC_COMM_WORLD, "# Mesh resolution: %f\n", opt->xc[1] / (opt->nxyz[0] - 1));
    PetscPrintf(PETSC_COMM_WORLD, "##############################################################\n");


    // STEP 7: CLEAN UP AFTER YOURSELF
    delete mma;
    delete output;
    delete filter;
    delete opt;
    delete physics;

    // Finalize PETSc / MPI
    PetscFinalize();

    return 1;
}
