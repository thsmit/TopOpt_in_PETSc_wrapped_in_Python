#include "topoptlib.h" // for wrapper
#include "Filter.h"
#include "LinearElasticity.h"
#include "MMA.h"
#include "MPIIO.h"
#include "TopOpt.h"
#include "mpi.h"
#include <petsc.h>

/*
Modified: Thijs Smit, May 2020

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

    // STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
    TopOpt* opt = new TopOpt(data);

    // STEP 2: THE PHYSICS
    LinearElasticity* physics = new LinearElasticity(opt->da_nodes, data);

    // STEP 3: THE FILTERING
    //Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->filter, opt->rmin);
    Filter* filter = new Filter(opt->da_nodes, opt->xPhys, opt->xPassive, opt->filter, opt->rmin);

    // STEP 4: VISUALIZATION USING VTK
    MPIIO* output = new MPIIO(opt->da_nodes, 3, "ux, uy, uz", 3, "x, xTilde, xPhys");
    
    // STEP 5: THE OPTIMIZER MMA
    MMA* mma;
    //MMA* mma = new MMA(data.nael, 1, opt->xMMA);
    PetscInt itr = 0;
    opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !
    // mma->SetAsymptotes(0.2, 0.65, 1.05);

    // STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
    ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
    CHKERRQ(ierr);

    // STEP 7: OPTIMIZATION LOOP
    PetscScalar ch = 1.0;
    double      t1, t2;
    while (itr < opt->maxItr && ch > 0.01) {
        // Update iteration counter
        itr++;

        // UPDATING THE CONTINUATION PARAMETERS
        // reference, Maximum Size...
		if (opt->continuationStatus && itr % opt->IterProj == 0){
			opt->penal   = PetscMin(opt->penal+0.25,3); // penal = penal : 0.25 : 3.0
			// move limits: initial 0.6, final 0.05
			//opt->movlim  = (opt->movlimEnd-opt->movlimIni)/(3.0-1.0)*(opt->penal-1.0)+opt->movlimIni; 
			//opt->Beta    = PetscMin((1.50*opt->Beta),38.0); // beta = beta*1.5 
			//PetscPrintf(PETSC_COMM_WORLD,"===================================================\n");
			PetscPrintf(PETSC_COMM_WORLD,"It.: %i, Penal: %2.2f   movlim: %f\n", itr, opt->penal, opt->movlim);
			//PetscPrintf(PETSC_COMM_WORLD,"===================================================\n");
		}

        // start timer
        t1 = MPI_Wtime();

        // Compute (a) obj+const, (b) sens, (c) obj+const+sens
        // input -> xPhys + SIMP settings + material properties
        // output -> objective + sensitivities + constraint value + cons sensitivies
        // output -> fx, dfdx, gx, dgdx
        ierr = physics->ComputeObjectiveConstraintsSensitivities(&(opt->fx), &(opt->gx[0]), opt->dfdx, opt->dgdx[0],
                                                                 opt->xPhys, opt->Emin, opt->Emax, opt->penal,
                                                                 opt->volfrac, data);

        
        CHKERRQ(ierr);

        // Compute objective scale
        if (itr == 1 || (itr % opt->IterProj == 0 && opt->continuationStatus)) {
            opt->fscale = 10.0 / opt->fx;
        }
        // Scale objectie and sens
        opt->fx = opt->fx * opt->fscale;
        VecScale(opt->dfdx, opt->fscale);

        // Filter sensitivities (chainrule)
        ierr = filter->Gradients(opt->x, opt->xTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, opt->beta,
                                 opt->eta);
        CHKERRQ(ierr);

        if (opt->xPassiveStatus) {
            
            //PetscPrintf(PETSC_COMM_WORLD, "just checking\n");
            
            // map vectors
            opt->UpdateVariables(1, opt->x, opt->xMMA);
            opt->UpdateVariables(1, opt->dfdx, opt->dfdxMMA);
            opt->UpdateVariables(1, opt->dgdx[0], opt->dgdxMMA[0]);
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
            opt->UpdateVariables(-1, opt->dgdx[0], opt->dgdxMMA[0]);
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
            changeBeta = filter->IncreaseBeta(&(opt->beta), opt->betaFinal, opt->gx[0], itr, ch);
        }

        // Filter design field
        ierr = filter->FilterProject(opt->x, opt->xTilde, opt->xPhys, opt->projectionFilter, opt->beta, opt->eta);
        CHKERRQ(ierr);

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
            output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr);
        }

        // Dump data needed for restarting code at termination
        if (itr % 10 == 0) {
            opt->WriteRestartFiles(&itr, mma);
            physics->WriteRestartFiles();
        }
    }
        
    // Write restart WriteRestartFiles
    opt->WriteRestartFiles(&itr, mma);
    physics->WriteRestartFiles();

    // Dump final design
    output->WriteVTK(physics->da_nodal, physics->GetStateField(), opt->x, opt->xTilde, opt->xPhys, itr + 1);

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
