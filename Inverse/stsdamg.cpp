 
#include "petscda.h"            /*I "petscda.h"   I*/
#include "petscksp.h"           /*I "petscksp.h"  I*/
#include "petscmg.h"            /*I "petscmg.h"   I*/
#include "stsdamgHeader.h"          
#include "private/pcimpl.h"  /*I "petscpc.h"   I*/

/*
   Code for almost fully managing multigrid/multi-level linear solvers for DA grids
*/

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGCreate"
/*@C
    stsDMMGCreate - Creates a DA based multigrid solver object. This allows one to 
      easily implement MG methods on regular grids.

    Collective on MPI_Comm

    Input Parameter:
+   comm - the processors that will share the grids and solution process
.   nlevels - number of multigrid levels 
-   user - an optional user context

    Output Parameters:
.    - the context

    Notes:
      To provide a different user context for each level call stsDMMGSetUser() after calling
      this routine

    Level: advanced

.seealso stsDMMGDestroy(), stsDMMGSetUser(), stsDMMGGetUser()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGCreate(MPI_Comm comm,PetscInt nlevels,void *user,stsDMMG **dmmg)
{
  PetscErrorCode ierr;
  PetscInt       i;
  stsDMMG           *p;
  PetscTruth     galerkin=PETSC_FALSE;

  PetscFunctionBegin;
  ierr = PetscOptionsGetInt(0,"-dmmg_nlevels",&nlevels,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(0,"-dmmg_galerkin",&galerkin);CHKERRQ(ierr);

  ierr = PetscMalloc(nlevels*sizeof(stsDMMG),&p);CHKERRQ(ierr);
  for (i=0; i<nlevels; i++) {
    ierr           = PetscNew(struct _p_stsDMMG,&p[i]);CHKERRQ(ierr);
    ierr           = PetscMemzero(p[i],sizeof(struct _p_stsDMMG));CHKERRQ(ierr);
    p[i]->nlevels  = nlevels - i;
    p[i]->comm     = comm;
    p[i]->user     = user;
    p[i]->galerkin = galerkin;
    //matFreeInterpolation Routines
    p[i]->createInterpolationMF = 0;
    p[i]->computeInterpolationMF = 0;
  }
  p[nlevels-1]->galerkin = PETSC_FALSE;
  *dmmg = p;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetUseGalerkinCoarse"
/*@C
    stsDMMGSetUseGalerkinCoarse - Courses the stsDMMG to use R*A_f*R^T to form
       the coarser matrices from finest 

    Collective on stsDMMG

    Input Parameter:
.    - the context

    Options Database Keys:
.    -dmmg_galerkin

    Level: advanced

    Notes: After you have called this you can manually set dmmg[0]->galerkin = PETSC_FALSE
       to have the coarsest grid not compute via Galerkin but still have the intermediate
       grids computed via Galerkin.

.seealso stsDMMGCreate()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUseGalerkinCoarse(stsDMMG* dmmg)
{
  PetscInt  i,nlevels = dmmg[0]->nlevels;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");

  for (i=0; i<nlevels-1; i++) {
    dmmg[i]->galerkin = PETSC_TRUE;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGDestroy"
/*@C
    stsDMMGDestroy - Destroys a DA based multigrid solver object. 

    Collective on stsDMMG

    Input Parameter:
.    - the context

    Level: advanced

.seealso stsDMMGCreate()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGDestroy(stsDMMG *dmmg)
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");

  for (i=1; i<nlevels; i++) {
    if (dmmg[i]->R) {ierr = MatDestroy(dmmg[i]->R);CHKERRQ(ierr);}
  }
  for (i=0; i<nlevels; i++) {
    if (dmmg[i]->dm)      {ierr = DMDestroy(dmmg[i]->dm);CHKERRQ(ierr);}
    if (dmmg[i]->x)       {ierr = VecDestroy(dmmg[i]->x);CHKERRQ(ierr);}
    if (dmmg[i]->b)       {ierr = VecDestroy(dmmg[i]->b);CHKERRQ(ierr);}
    if (dmmg[i]->r)       {ierr = VecDestroy(dmmg[i]->r);CHKERRQ(ierr);}
    if (dmmg[i]->work1)   {ierr = VecDestroy(dmmg[i]->work1);CHKERRQ(ierr);}
    if (dmmg[i]->w)       {ierr = VecDestroy(dmmg[i]->w);CHKERRQ(ierr);}
    if (dmmg[i]->work2)   {ierr = VecDestroy(dmmg[i]->work2);CHKERRQ(ierr);}
    if (dmmg[i]->lwork1)  {ierr = VecDestroy(dmmg[i]->lwork1);CHKERRQ(ierr);}
    if (dmmg[i]->B && dmmg[i]->B != dmmg[i]->J) {ierr = MatDestroy(dmmg[i]->B);CHKERRQ(ierr);}
    if (dmmg[i]->J)         {ierr = MatDestroy(dmmg[i]->J);CHKERRQ(ierr);}
    if (dmmg[i]->Rscale)    {ierr = VecDestroy(dmmg[i]->Rscale);CHKERRQ(ierr);}
    if (dmmg[i]->fdcoloring){ierr = MatFDColoringDestroy(dmmg[i]->fdcoloring);CHKERRQ(ierr);}
    if (dmmg[i]->ksp && !dmmg[i]->snes) {ierr = KSPDestroy(dmmg[i]->ksp);CHKERRQ(ierr);}
    if (dmmg[i]->snes)      {ierr = PetscObjectDestroy((PetscObject)dmmg[i]->snes);CHKERRQ(ierr);} 
    if (dmmg[i]->inject)    {ierr = VecScatterDestroy(dmmg[i]->inject);CHKERRQ(ierr);} 
    ierr = PetscFree(dmmg[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dmmg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetDM"
/*@C
    stsDMMGSetDM - Sets the coarse grid information for the grids

    Collective on stsDMMG

    Input Parameter:
+   dmmg - the context
-   dm - the DA or VecPack object

    Level: advanced

.seealso stsDMMGCreate(), stsDMMGDestroy()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetDM(stsDMMG *dmmg,DM dm)
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");

  /* Create DA data structure for all the levels */
  dmmg[0]->dm = dm;
  ierr = PetscObjectReference((PetscObject)dm);CHKERRQ(ierr);
  for (i=1; i<nlevels; i++) {
    ierr = DMRefine(dmmg[i-1]->dm,dmmg[i]->comm,&dmmg[i]->dm);CHKERRQ(ierr);
  }
  ierr = stsDMMGSetUp(dmmg);CHKERRQ(ierr); 
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetUp"
/*@C
    stsDMMGSetUp - Prepares the stsDMMG to solve a system

    Collective on stsDMMG

    Input Parameter:
.   dmmg - the context

    Level: advanced

.seealso stsDMMGCreate(), stsDMMGDestroy(), stsDMMG, stsDMMGSetSNES(), stsDMMGSetKSP(), stsDMMGSolve()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUp(stsDMMG *dmmg)
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;
  PetscInt resType = 0;//R= P' using Petsc's stardard operators
  
  PetscFunctionBegin;

  /* Create work vectors and matrix for each level */
  for (i=0; i<nlevels; i++) {
    ierr = DMCreateGlobalVector(dmmg[i]->dm,&dmmg[i]->x); CHKERRQ(ierr);
    ierr = VecDuplicate(dmmg[i]->x,&dmmg[i]->b); CHKERRQ(ierr);
    ierr = VecDuplicate(dmmg[i]->x,&dmmg[i]->r); CHKERRQ(ierr);
  }

  
  ierr = PetscOptionsGetInt(0,"-restype",&resType,0); CHKERRQ(ierr);
  /* Create interpolation/restriction between levels.\
     resType =0 for std. petsc routines, resType .ne. 0 for Mat-Free routines */
  if(resType) {
    if((nlevels>1) && (!dmmg[0]->createInterpolationMF || !dmmg[0]->computeInterpolationMF)) {
      SETERRQ(PETSC_ERR_USER,"Routines for creating and computing Matrix-Free Interpolations are not provided.");
    }
    for (i=1; i<nlevels; i++) {
      //nlevels must be atleast 2 for Restriction/Interpolation
      /*This has been changed from DMGetInterpolation inorder to use matrix-free Interpolation.
	Note: At present, the matrix-free equivalent method has been implemented only for 3-D structured grids with
	1-dof/node and nonperiodic B.C. For all other cases, the standard DA routines will be used. Also, Note.
	It has been assumed that the DM object is of type DA and not VecPac. */
      ierr = (*dmmg[i-1]->createInterpolationMF)(dmmg[i-1],dmmg[i],&dmmg[i]->R); CHKERRQ(ierr);
      ierr = (*dmmg[i-1]->computeInterpolationMF)(dmmg[i-1],dmmg[i],dmmg[i]->R ); CHKERRQ(ierr);
    }//end for
  }else {
    for (i=1; i<nlevels; i++) {
      printf("Inside CreateInterpolation using Std. Petsc. Routines. Level = %d\n",dmmg[i-1]->nlevels);
      ierr = DMGetInterpolation(dmmg[i-1]->dm,dmmg[i]->dm,&dmmg[i]->R,PETSC_NULL); CHKERRQ(ierr);//RS:Note no scaling
      printf("Finished CreateInterpolation. using Std. Petsc. Routines. Level = %d \n",dmmg[i-1]->nlevels);
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "stsDMMGSetInterpolationMatrixFree"
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetInterpolationMatrixFree(stsDMMG *dmmg,
     PetscErrorCode (*crIntrp)(stsDMMG,stsDMMG,Mat*),PetscErrorCode (*compIntrp)(stsDMMG,stsDMMG,Mat)) {
    
  PetscInt       i,nlevels = dmmg[0]->nlevels;
  PetscFunctionBegin;
  for (i=0; i<nlevels-1; i++) {
    dmmg[i]->createInterpolationMF = crIntrp;
    dmmg[i]->computeInterpolationMF = compIntrp;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "stsDMMGCreateJMatrix"
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGCreateJMatrix(stsDMMG dmmg,PetscErrorCode (*crjac)(stsDMMG,Mat*)) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(crjac){
        ierr = (*crjac)(dmmg,&(dmmg->J)); CHKERRQ(ierr);
  }else{
        ierr = DMGetMatrix(dmmg->dm,MATAIJ,&dmmg->J); CHKERRQ(ierr);//RS:This calls MatCreate!
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSolve"
/*@C
    stsDMMGSolve - Actually solves the (non)linear system defined with the stsDMMG

    Collective on stsDMMG

    Input Parameter:
.   dmmg - the context

    Level: advanced

    Options Database:
+   -dmmg_grid_sequence - use grid sequencing to get the initial solution for each level from the previous
-   -dmmg_vecmonitor - display the solution at each iteration

     Notes: For linear (KSP) problems may be called more than once, uses the same 
    matrices but recomputes the right hand side for each new solve. Call stsDMMGSetKSP()
    to generate new matrices.
 
.seealso stsDMMGCreate(), stsDMMGDestroy(), stsDMMG, stsDMMGSetSNES(), stsDMMGSetKSP(), stsDMMGSetUp()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSolve(stsDMMG *dmmg)
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;
  PetscTruth     gridseq,vecmonitor,flg;

  PetscFunctionBegin;
  ierr = PetscOptionsHasName(0,"-dmmg_grid_sequence",&gridseq);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(0,"-dmmg_vecmonitor",&vecmonitor);CHKERRQ(ierr);
  if (gridseq) {
    if (dmmg[0]->initialguess) {
      ierr = (*dmmg[0]->initialguess)(dmmg[0],dmmg[0]->x);CHKERRQ(ierr);
      if (dmmg[0]->ksp && !dmmg[0]->snes) {
        ierr = KSPSetInitialGuessNonzero(dmmg[0]->ksp,PETSC_TRUE);CHKERRQ(ierr);
      }
    }
    for (i=0; i<nlevels-1; i++) {
      ierr = (*dmmg[i]->solve)(dmmg,i);CHKERRQ(ierr);
      if (vecmonitor) {
        ierr = VecView(dmmg[i]->x,PETSC_VIEWER_DRAW_(dmmg[i]->comm));CHKERRQ(ierr);
      }
      ierr = MatInterpolate(dmmg[i+1]->R,dmmg[i]->x,dmmg[i+1]->x);CHKERRQ(ierr);
      if (dmmg[i+1]->ksp && !dmmg[i+1]->ksp) {
        ierr = KSPSetInitialGuessNonzero(dmmg[i+1]->ksp,PETSC_TRUE);CHKERRQ(ierr);
      }
    }
  } else {
    if (dmmg[nlevels-1]->initialguess) {
      ierr = (*dmmg[nlevels-1]->initialguess)(dmmg[nlevels-1],dmmg[nlevels-1]->x);CHKERRQ(ierr);
    }
  }
  ierr = (*stsDMMGGetFine(dmmg)->solve)(dmmg,nlevels-1);CHKERRQ(ierr);
  if (vecmonitor) {
     ierr = VecView(dmmg[nlevels-1]->x,PETSC_VIEWER_DRAW_(dmmg[nlevels-1]->comm));CHKERRQ(ierr);
  }

  ierr = PetscOptionsHasName(PETSC_NULL,"-dmmg_view",&flg);CHKERRQ(ierr);
  if (flg && !PetscPreLoadingOn) {
    ierr = stsDMMGView(dmmg,PETSC_VIEWER_STDOUT_(dmmg[0]->comm));CHKERRQ(ierr);
  }
  ierr = PetscOptionsHasName(PETSC_NULL,"-dmmg_view_binary",&flg);CHKERRQ(ierr);
  if (flg && !PetscPreLoadingOn) {
    ierr = stsDMMGView(dmmg,PETSC_VIEWER_BINARY_(dmmg[0]->comm));CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSolveKSP"
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSolveKSP(stsDMMG *dmmg,PetscInt level)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (dmmg[level]->rhs) {
    ierr = (*dmmg[level]->rhs)(dmmg[level],dmmg[level]->b);CHKERRQ(ierr); 
  }
  if (dmmg[level]->matricesset) {
    ierr = KSPSetOperators(dmmg[level]->ksp,dmmg[level]->J,dmmg[level]->B,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
    dmmg[level]->matricesset = PETSC_FALSE;
  }
  ierr = KSPSolve(dmmg[level]->ksp,dmmg[level]->b,dmmg[level]->x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
    Sets each of the linear solvers to use multigrid 
*/
#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetUpLevel"
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUpLevel(stsDMMG *dmmg,KSP ksp,PetscInt nlevels)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PC             pc;
  PetscTruth     ismg,monitor,ismf,isshell,ismffd;
  KSP            lksp; /* solver internal to the multigrid preconditioner */
  MPI_Comm       *comms,comm;
  PetscViewer    ascii;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");

  ierr = PetscOptionsHasName(PETSC_NULL,"-dmmg_ksp_monitor",&monitor);CHKERRQ(ierr);
  if (monitor) {
    ierr = PetscObjectGetComm((PetscObject)ksp,&comm);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(comm,"stdout",&ascii);CHKERRQ(ierr);
    ierr = PetscViewerASCIISetTab(ascii,1+dmmg[0]->nlevels-nlevels);CHKERRQ(ierr);
    ierr = KSPMonitorSet(ksp,KSPMonitorDefault,ascii,(PetscErrorCode(*)(void*))PetscViewerDestroy);CHKERRQ(ierr);
  }

  /* use fgmres on outer iteration by default */
  ierr  = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
  ierr  = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr  = PCSetType(pc,PCMG);CHKERRQ(ierr);
  ierr  = PetscMalloc(nlevels*sizeof(MPI_Comm),&comms);CHKERRQ(ierr);
  for (i=0; i<nlevels; i++) {
    comms[i] = dmmg[i]->comm;
  }
  ierr  = PCMGSetLevels(pc,nlevels,comms);CHKERRQ(ierr);
  ierr  = PetscFree(comms);CHKERRQ(ierr); 
  ierr =  PCMGSetType(pc,PC_MG_FULL);CHKERRQ(ierr);

  ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
  if (ismg) {
    if (dmmg[0]->galerkin) {
      ierr = PCMGSetGalerkin(pc);CHKERRQ(ierr);
    }

    /* set solvers for each level */
    for (i=0; i<nlevels; i++) {
      ierr = PCMGGetSmoother(pc,i,&lksp);CHKERRQ(ierr);
      if (i == nlevels-1 || !dmmg[0]->galerkin) {
        ierr = KSPSetOperators(lksp,dmmg[i]->J,dmmg[i]->B,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
      }
      if (i < nlevels-1) { /* don't set for finest level, they are set in PCApply_MG()*/
	ierr = PCMGSetX(pc,i,dmmg[i]->x);CHKERRQ(ierr); 
	ierr = PCMGSetRhs(pc,i,dmmg[i]->b);CHKERRQ(ierr); 
      }
      if (i > 0) {
        ierr = PCMGSetR(pc,i,dmmg[i]->r);CHKERRQ(ierr); 
        ierr = PCMGSetResidual(pc,i,PCMGDefaultResidual,dmmg[i]->J);CHKERRQ(ierr);
      }
      if (monitor) {
        ierr = PetscObjectGetComm((PetscObject)lksp,&comm);CHKERRQ(ierr);
        ierr = PetscViewerASCIIOpen(comm,"stdout",&ascii);CHKERRQ(ierr);
        ierr = PetscViewerASCIISetTab(ascii,1+dmmg[0]->nlevels-i);CHKERRQ(ierr);
        ierr = KSPMonitorSet(lksp,KSPMonitorDefault,ascii,(PetscErrorCode(*)(void*))PetscViewerDestroy);CHKERRQ(ierr);
      }
      /* If using a matrix free multiply and did not provide an explicit matrix to build
         the preconditioner then must use no preconditioner 
      */
      ierr = PetscTypeCompare((PetscObject)dmmg[i]->B,MATSHELL,&isshell);CHKERRQ(ierr);
      ierr = PetscTypeCompare((PetscObject)dmmg[i]->B,MATDAAD,&ismf);CHKERRQ(ierr);
      ierr = PetscTypeCompare((PetscObject)dmmg[i]->B,MATMFFD,&ismffd);CHKERRQ(ierr);
      if (isshell || ismf || ismffd) {
        PC  lpc;
        ierr = KSPGetPC(lksp,&lpc);CHKERRQ(ierr);
        ierr = PCSetType(lpc,PCNONE);CHKERRQ(ierr);
      }
    }

    /* Set interpolation/restriction between levels */
    for (i=1; i<nlevels; i++) {
      ierr = PCMGSetInterpolation(pc,i,dmmg[i]->R);CHKERRQ(ierr); 
      ierr = PCMGSetRestriction(pc,i,dmmg[i]->R);CHKERRQ(ierr); 
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetKSP"
/*@C
    stsDMMGSetKSP - Sets the linear solver object that will use the grid hierarchy

    Collective on stsDMMG

    Input Parameter:
+   dmmg - the context
.   func - function to compute linear system matrix on each grid level
-   rhs - function to compute right hand side on each level (need only work on the finest grid
          if you do not use grid sequencing)

    Level: advanced

    Notes: For linear problems my be called more than once, reevaluates the matrices if it is called more
       than once. Call stsDMMGSolve() directly several times to solve with the same matrix but different 
       right hand sides.
   
.seealso stsDMMGCreate(), stsDMMGDestroy, stsDMMGSetDM(), stsDMMGSolve()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetKSP(stsDMMG *dmmg, PetscErrorCode (*crjac)(stsDMMG,Mat*),
	PetscErrorCode (*rhs)(stsDMMG,Vec),PetscErrorCode (*func)(stsDMMG,Mat,Mat),
	PetscErrorCode (*coarsefunc)(stsDMMG,Mat,Mat))
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;
  PetscTruth     galerkin;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");
  if(crjac != PETSC_NULL) {
    //ignore galerkin for matrix-free method.
    for(i=0;i<nlevels;i++) {
      dmmg[i]->galerkin = PETSC_FALSE;
    }
  }
  galerkin = dmmg[nlevels - 2 > 0 ? nlevels - 2 : 0]->galerkin;  

  if (galerkin) {
    ierr = DMGetMatrix(dmmg[nlevels-1]->dm,MATAIJ,&dmmg[nlevels-1]->B);CHKERRQ(ierr);
    if (!dmmg[nlevels-1]->J) {
      dmmg[nlevels-1]->J = dmmg[nlevels-1]->B;
    }
    ierr = (*func)(dmmg[nlevels-1],dmmg[nlevels-1]->J,dmmg[nlevels-1]->B);CHKERRQ(ierr);
    for (i=nlevels-2; i>-1; i--) {
      if (dmmg[i]->galerkin) {
        ierr = MatPtAP(dmmg[i+1]->B,dmmg[i+1]->R,MAT_INITIAL_MATRIX,1.0,&dmmg[i]->B);CHKERRQ(ierr);
        if (!dmmg[i]->J) {
          dmmg[i]->J = dmmg[i]->B;
        }
      }
    }
  }

  if (!dmmg[0]->ksp) {
    /* create solvers for each level if they don't already exist*/
    for (i=0; i<nlevels; i++) {

      if (!dmmg[i]->B && !dmmg[i]->galerkin) {
	// GB: this is were the matrix is created!
	if( (crjac == PETSC_NULL) || ( (i==0) && (coarsefunc != PETSC_NULL) ) ){
	  ierr = DMGetMatrix(dmmg[i]->dm,MATAIJ,&dmmg[i]->B);CHKERRQ(ierr);//RS:This calls MatCreate!
	  //So for the Matrix-free method, crjac must not be NULL!
	  //That means we must provide a create Jacobian method!
	}
	else{
	  ierr = (*crjac)(dmmg[i],&dmmg[i]->B);CHKERRQ(ierr);
	}//end if-else

      } 
      if (!dmmg[i]->J) {
        dmmg[i]->J = dmmg[i]->B;
      }

      ierr = KSPCreate(dmmg[i]->comm,&dmmg[i]->ksp);CHKERRQ(ierr);
      ierr = stsDMMGSetUpLevel(dmmg,dmmg[i]->ksp,i+1);CHKERRQ(ierr);
      ierr = KSPSetFromOptions(dmmg[i]->ksp);CHKERRQ(ierr);
      dmmg[i]->solve = stsDMMGSolveKSP;
      dmmg[i]->rhs   = rhs;
    }
  }

  /* evalute matrix on each level */
 if(!dmmg[0]->galerkin){
    if(coarsefunc == PETSC_NULL) {
      ierr = (*func)(dmmg[0],dmmg[0]->J,dmmg[0]->B);CHKERRQ(ierr);
    }else {
      ierr = (*coarsefunc)(dmmg[0],dmmg[0]->J,dmmg[0]->B);CHKERRQ(ierr);
    }//end if-else
 }  
    dmmg[0]->matricesset = PETSC_TRUE;
  for (i=1; i<nlevels; i++) {
    if (!dmmg[i]->galerkin) {
      ierr = (*func)(dmmg[i],dmmg[i]->J,dmmg[i]->B);CHKERRQ(ierr);
    }
    dmmg[i]->matricesset = PETSC_TRUE;
  }

  for (i=0; i<nlevels-1; i++) {
    ierr = KSPSetOptionsPrefix(dmmg[i]->ksp,"dmmg_");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGView"
/*@C
    stsDMMGView - prints information on a DA based multi-level preconditioner

    Collective on stsDMMG and PetscViewer

    Input Parameter:
+   dmmg - the context
-   viewer - the viewer

    Level: advanced

.seealso stsDMMGCreate(), stsDMMGDestroy

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGView(stsDMMG *dmmg,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       i,nlevels = dmmg[0]->nlevels;
  PetscMPIInt    flag;
  MPI_Comm       comm;
  PetscTruth     iascii,isbinary;

  PetscFunctionBegin;
  PetscValidPointer(dmmg,1);
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_COOKIE,2);
  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_compare(comm,dmmg[0]->comm,&flag);CHKERRQ(ierr);
  if (flag != MPI_CONGRUENT && flag != MPI_IDENT) {
    SETERRQ(PETSC_ERR_ARG_NOTSAMECOMM,"Different communicators in the stsDMMG and the PetscViewer");
  }

  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii);CHKERRQ(ierr);
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_BINARY,&isbinary);CHKERRQ(ierr);
  if (isbinary) {
    for (i=0; i<nlevels; i++) {
      ierr = MatView(dmmg[i]->J,viewer);CHKERRQ(ierr);
    }
    for (i=1; i<nlevels; i++) {
      ierr = MatView(dmmg[i]->R,viewer);CHKERRQ(ierr);
    }
  } else {
    if (iascii) {
      ierr = PetscViewerASCIIPrintf(viewer,"stsDMMG Object with %D levels\n",nlevels);CHKERRQ(ierr);
    }
    for (i=0; i<nlevels; i++) {
      ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
      ierr = DMView(dmmg[i]->dm,viewer);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
    }
    if (iascii) {
      ierr = PetscViewerASCIIPrintf(viewer,"%s Object on finest level\n",dmmg[nlevels-1]->ksp ? "KSP" : "SNES");CHKERRQ(ierr);
      if (dmmg[nlevels-2 > 0 ? nlevels-2 : 0]->galerkin) {
	ierr = PetscViewerASCIIPrintf(viewer,"Using Galerkin R^T*A*R process to compute coarser matrices");CHKERRQ(ierr);
      }
    }
    if (dmmg[nlevels-1]->ksp) {
      ierr = KSPView(dmmg[nlevels-1]->ksp,viewer);CHKERRQ(ierr);
    } else {
      /* use of PetscObjectView() means we do not have to link with libpetscsnes if SNES is not being used */
      ierr = PetscObjectView((PetscObject)dmmg[nlevels-1]->snes,viewer);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetNullSpace"
/*@C
    stsDMMGSetNullSpace - Indicates the null space in the linear operator (this is needed by the linear solver)

    Collective on stsDMMG

    Input Parameter:
+   dmmg - the context
.   has_cnst - is the constant vector in the null space
.   n - number of null vectors (excluding the possible constant vector)
-   func - a function that fills an array of vectors with the null vectors (must be orthonormal), may be PETSC_NULL

    Level: advanced

.seealso stsDMMGCreate(), stsDMMGDestroy, stsDMMGSetDM(), stsDMMGSolve(), MatNullSpaceCreate(), KSPSetNullSpace()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetNullSpace(stsDMMG *dmmg,PetscTruth has_cnst,PetscInt n,PetscErrorCode (*func)(stsDMMG,Vec[]))
{
  PetscErrorCode ierr;
  PetscInt       i,j,nlevels = dmmg[0]->nlevels;
  Vec            *nulls = 0;
  MatNullSpace   nullsp;
  KSP            iksp;
  PC             pc,ipc;
  PetscTruth     ismg,isred;

  PetscFunctionBegin;
  if (!dmmg) SETERRQ(PETSC_ERR_ARG_NULL,"Passing null as stsDMMG");
  if (!dmmg[0]->ksp) SETERRQ(PETSC_ERR_ORDER,"Must call AFTER stsDMMGSetKSP() or stsDMMGSetSNES()");
  if ((n && !func) || (!n && func)) SETERRQ(PETSC_ERR_ARG_INCOMP,"Both n and func() must be set together");
  if (n < 0) SETERRQ1(PETSC_ERR_ARG_OUTOFRANGE,"Cannot have negative number of vectors in null space n = %D",n)

  for (i=0; i<nlevels; i++) {
    if (n) {
      ierr = VecDuplicateVecs(dmmg[i]->b,n,&nulls);CHKERRQ(ierr);
      ierr = (*func)(dmmg[i],nulls);CHKERRQ(ierr);
    }
    ierr = MatNullSpaceCreate(dmmg[i]->comm,has_cnst,n,nulls,&nullsp);CHKERRQ(ierr);
    ierr = KSPSetNullSpace(dmmg[i]->ksp,nullsp);CHKERRQ(ierr);
    for (j=i; j<nlevels; j++) {
      ierr = KSPGetPC(dmmg[j]->ksp,&pc);CHKERRQ(ierr);
      ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
      if (ismg) {
        ierr = PCMGGetSmoother(pc,i,&iksp);CHKERRQ(ierr);
        ierr = KSPSetNullSpace(iksp, nullsp);CHKERRQ(ierr);
      }
    }
    ierr = MatNullSpaceDestroy(nullsp);CHKERRQ(ierr);
    if (n) {
      ierr = PetscFree(nulls);CHKERRQ(ierr);
    }
  }
  /* make all the coarse grid solvers have LU shift since they are singular */
  for (i=0; i<nlevels; i++) {
    ierr = KSPGetPC(dmmg[i]->ksp,&pc);CHKERRQ(ierr);
    ierr = PetscTypeCompare((PetscObject)pc,PCMG,&ismg);CHKERRQ(ierr);
    if (ismg) {
      ierr = PCMGGetSmoother(pc,0,&iksp);CHKERRQ(ierr);
      ierr = KSPGetPC(iksp,&ipc);CHKERRQ(ierr);
      ierr = PetscTypeCompare((PetscObject)ipc,PCREDUNDANT,&isred);CHKERRQ(ierr);
      if (isred) {
        ierr = PCRedundantGetPC(ipc,&ipc);CHKERRQ(ierr);
      }
      ierr = PCFactorSetShiftPd(ipc,PETSC_TRUE);CHKERRQ(ierr); 
    }
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGInitialGuessCurrent"
/*@C
    stsDMMGInitialGuessCurrent - Use with stsDMMGSetInitialGuess() to use the current value in the 
       solution vector (obtainable with stsDMMGGetx() as the initial guess)

    Collective on stsDMMG

    Input Parameter:
+   dmmg - the context
-   vec - dummy argument

    Level: intermediate

.seealso stsDMMGCreate(), stsDMMGDestroy, stsDMMGSetKSP(), stsDMMGSetSNES(), stsDMMGSetInitialGuess()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGInitialGuessCurrent(stsDMMG dmmg,Vec vec)
{
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "stsDMMGSetInitialGuess"
/*@C
    stsDMMGSetInitialGuess - Sets the function that computes an initial guess.

    Collective on stsDMMG

    Input Parameter:
+   dmmg - the context
-   guess - the function

    Notes: For nonlinear problems, if this is not set, then the current value in the 
             solution vector (obtained with stsDMMGGetX()) is used. Thus is if you doing 'time
             stepping' it will use your current solution as the guess for the next timestep.
           If grid sequencing is used (via -dmmg_grid_sequence) then the "guess" function
             is used only on the coarsest grid.
           For linear problems, if this is not set, then 0 is used as an initial guess.
             If you would like the linear solver to also (like the nonlinear solver) use
             the current solution vector as the initial guess then use stsDMMGInitialGuessCurrent()
             as the function you pass in

    Level: intermediate


.seealso stsDMMGCreate(), stsDMMGDestroy, stsDMMGSetKSP(), stsDMMGSetSNES(), stsDMMGInitialGuessCurrent()

@*/
PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetInitialGuess(stsDMMG *dmmg,PetscErrorCode (*guess)(stsDMMG,Vec))
{
  PetscInt i,nlevels = dmmg[0]->nlevels;

  PetscFunctionBegin;
  for (i=0; i<nlevels; i++) {
    dmmg[i]->initialguess = guess;
  }
  PetscFunctionReturn(0);
}






