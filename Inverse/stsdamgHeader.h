/*
  Defines the interface functions for the stsDMMG object.
*/
#ifndef __STSDMMGHEADER_H
#define __STSDMMGHEADER_H
#include "petscsnes.h"
#include "petscda.h"
PETSC_EXTERN_CXX_BEGIN

/*S
     stsDMMG -  Data structure to easily manage multi-level non-linear solvers on grids managed by DM
          
   Level: intermediate

  Concepts: multigrid, Newton-multigrid

.seealso:  VecPackCreate(), DA, VecPack, DM, stsDMMGCreate(), stsDMMGSetKSP(), stsDMMGSetSNES()
S*/
typedef struct _p_stsDMMG* stsDMMG;
struct _p_stsDMMG {
  DM             dm;                   /* grid information for this level */
  Vec            x,b,r;                /* global vectors used in multigrid preconditioner for this level*/
  Mat            J;                    /* matrix on this level */
  Mat            B;
  Mat            R;                    /* interpolation to next finer level (actually a misnomer) */
  PetscInt       nlevels;              /* number of levels above this one (total number of levels on level 0)*/
  MPI_Comm       comm;
  PetscErrorCode (*solve)(stsDMMG*,PetscInt);
  void           *user;         
  PetscTruth     galerkin;                  /* for A_c = R*A*R^T */

  /*Matrix-Free Intergrid Transfer Operators */
  PetscErrorCode (*createInterpolationMF)(stsDMMG,stsDMMG,Mat*);
  PetscErrorCode (*computeInterpolationMF)(stsDMMG,stsDMMG,Mat);
  
  /* KSP only */
  KSP            ksp;             
  PetscErrorCode (*rhs)(stsDMMG,Vec);
  PetscTruth     matricesset;               /* User had called stsDMMGSetKSP() and the matrices have been computed */

  /* SNES only */
  Vec            Rscale;                 /* scaling to restriction before computing Jacobian */
  PetscErrorCode (*computejacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*);  
  PetscErrorCode (*computefunction)(SNES,Vec,Vec,void*);  

  PetscTruth     updatejacobian;         /* compute new Jacobian when stsDMMGComputeJacobian_Multigrid() is called */
  PetscInt       updatejacobianperiod;   /* how often, inside a SNES, the Jacobian is recomputed */

  MatFDColoring  fdcoloring;             /* only used with FD coloring for Jacobian */  
  SNES           snes;                  
  PetscErrorCode (*initialguess)(stsDMMG,Vec);
  Vec            w,work1,work2;         /* global vectors */
  Vec            lwork1;

  /* FAS only */
  NLF            nlf;                   /* FAS smoother object */
  VecScatter     inject;                /* inject from this level to the next coarsest */
  PetscTruth     monitor,monitorall;
  PetscInt       presmooth,postsmooth,coarsesmooth;
  PetscReal      rtol,abstol,rrtol;       /* convergence tolerance */   
  
};

//RS: This is the function used to set create and compute InterpolationMF. This must be called before calling stsDMMGSetUp
//Since stsDMMGSetUp is called in stsDMMGSetDM, this must be called before calling stsDMMGSetDM
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetInterpolationMatrixFree(stsDMMG*,
PetscErrorCode (*)(stsDMMG,stsDMMG,Mat*),PetscErrorCode (*)(stsDMMG,stsDMMG,Mat) );

//RS: This is the function used to create the J matrix. This must be called before calling stsDMMGSetKSP or stsDMMGSetSNES if you want to set
//J and B to be different. If this is not called, by default J and B will be equal. If the second argument is zero, a full-matrix will be
// created using DMGetMatrix(), else the user-provided function will be used to create the matrix (possibly a MatShell).
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGCreateJMatrix(stsDMMG,PetscErrorCode (*)(stsDMMG,Mat*));

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGCreate(MPI_Comm,PetscInt,void*,stsDMMG**);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGDestroy(stsDMMG*);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUp(stsDMMG*);

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetKSP(stsDMMG*,PetscErrorCode (*)(stsDMMG,Mat*),PetscErrorCode (*)(stsDMMG,Vec), PetscErrorCode (*)(stsDMMG,Mat,Mat), PetscErrorCode (*)(stsDMMG,Mat,Mat));

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetSNES(stsDMMG*,PetscErrorCode (*)(stsDMMG,Mat*),PetscErrorCode (*)(SNES,Vec,Vec,void*),PetscErrorCode (*)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),PetscErrorCode (*)(SNES,Vec,Mat*,Mat*,MatStructure*,void*));

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetInitialGuess(stsDMMG*,PetscErrorCode (*)(stsDMMG,Vec));
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGInitialGuessCurrent(stsDMMG,Vec);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGView(stsDMMG*,PetscViewer);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSolve(stsDMMG*);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUseMatrixFree(stsDMMG*);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetDM(stsDMMG*,DM);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUpLevel(stsDMMG*,KSP,PetscInt);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetUseGalerkinCoarse(stsDMMG*);
EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetNullSpace(stsDMMG*,PetscTruth,PetscInt,PetscErrorCode (*)(stsDMMG,Vec[]));

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetSNESLocal_Private(stsDMMG*,PetscErrorCode (*)(stsDMMG,Mat*),DALocalFunction1,DALocalFunction1,DALocalFunction1,DALocalFunction1,DALocalFunction1,DALocalFunction1);
#if defined(PETSC_HAVE_ADIC)
#  define stsDMMGSetSNESLocal(dmmg,crjac,function,jacobian,coarsejacobian,ad_function,ad_coarsefunction,admf_function) \
  stsDMMGSetSNESLocal_Private(dmmg,(PetscErrorCode (*)(stsDMMG,Mat*))crjac,\
  (DALocalFunction1)function,(DALocalFunction1)jacobian,(DALocalFunction1)coarsejacobian,\
  (DALocalFunction1)(ad_function),(DALocalFunction1)(ad_coarsefunction),(DALocalFunction1)(admf_function))
#else
#  define stsDMMGSetSNESLocal(dmmg,crjac,function,jacobian,coarsejacobian,ad_function,ad_coarsefunction,admf_function) \
    stsDMMGSetSNESLocal_Private(dmmg,(PetscErrorCode (*)(stsDMMG,Mat*))crjac,\
   (DALocalFunction1)function,(DALocalFunction1)jacobian,(DALocalFunction1)coarsejacobian,\
   (DALocalFunction1)0,(DALocalFunction1)0,(DALocalFunction1)0)
#endif

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetSNESLocali_Private(stsDMMG*,PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*),PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*),PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*));
#if defined(PETSC_HAVE_ADIC)
#  define stsDMMGSetSNESLocali(dmmg,function,ad_function,admf_function) stsDMMGSetSNESLocali_Private(dmmg,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*))function,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*))(ad_function),(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*))(admf_function))
#else
#  define stsDMMGSetSNESLocali(dmmg,function,ad_function,admf_function) stsDMMGSetSNESLocali_Private(dmmg,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*))function,0,0)
#endif

EXTERN PetscErrorCode PETSCSNES_DLLEXPORT stsDMMGSetSNESLocalib_Private(stsDMMG*,PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*),PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*),PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*));
#if defined(PETSC_HAVE_ADIC)
#  define stsDMMGSetSNESLocalib(dmmg,function,ad_function,admf_function) stsDMMGSetSNESLocalib_Private(dmmg,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*))function,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*))(ad_function),(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,void*,void*))(admf_function))
#else
#  define stsDMMGSetSNESLocalib(dmmg,function,ad_function,admf_function) stsDMMGSetSNESLocalib_Private(dmmg,(PetscErrorCode(*)(DALocalInfo*,MatStencil*,void*,PetscScalar*,void*))function,0,0)
#endif

/*MC
   stsDMMGGetRHS - Returns the right hand side vector from a stsDMMG solve on the finest grid

   Synopsis:
   Vec stsDMMGGetRHS(stsDMMG *dmmg)

   Not Collective, but resulting vector is parallel

   Input Parameters:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Fortran Usage:
.     stsDMMGGetRHS(stsDMMG dmmg,Vec b,PetscErrorCode ierr)

.seealso: stsDMMGCreate(), stsDMMGSetSNES(), stsDMMGSetKSP(), stsDMMGSetSNESLocal(), stsDMMGGetRHS()

M*/
#define stsDMMGGetRHS(ctx)              (ctx)[(ctx)[0]->nlevels-1]->b

#define stsDMMGGetr(ctx)              (ctx)[(ctx)[0]->nlevels-1]->r

/*MC
   stsDMMGGetx - Returns the solution vector from a stsDMMG solve on the finest grid

   Synopsis:
   Vec stsDMMGGetx(stsDMMG *dmmg)

   Not Collective, but resulting vector is parallel

   Input Parameters:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Fortran Usage:
.     stsDMMGGetx(stsDMMG dmmg,Vec x,PetscErrorCode ierr)

.seealso: stsDMMGCreate(), stsDMMGSetSNES(), stsDMMGSetKSP(), stsDMMGSetSNESLocal()

M*/
#define stsDMMGGetx(ctx)              (ctx)[(ctx)[0]->nlevels-1]->x

/*MC
   stsDMMGGetJ - Returns the Jacobian (matrix) for the finest level

   Synopsis:
   Mat stsDMMGGetJ(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetB(), stsDMMGGetRHS()

M*/
#define stsDMMGGetJ(ctx)              (ctx)[(ctx)[0]->nlevels-1]->J

/*MC
   stsDMMGGetComm - Returns the MPI_Comm for the finest level

   Synopsis:
   MPI_Comm stsDMMGGetJ(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ()

M*/
#define stsDMMGGetComm(ctx)           (ctx)[(ctx)[0]->nlevels-1]->comm

/*MC
   stsDMMGGetB - Returns the matrix for the finest level used to construct the preconditioner; usually 
              the same as the Jacobian

   Synopsis:
   Mat stsDMMGGetJ(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ()

M*/
#define stsDMMGGetB(ctx)              (ctx)[(ctx)[0]->nlevels-1]->B

/*MC
   stsDMMGGetFine - Returns the stsDMMG associated with the finest level

   Synopsis:
   stsDMMG stsDMMGGetFine(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ()

M*/
#define stsDMMGGetFine(ctx)           (ctx)[(ctx)[0]->nlevels-1]


/*MC
   stsDMMGGetKSP - Gets the KSP object (linear solver object) for the finest level

   Synopsis:
   KSP stsDMMGGetKSP(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Notes: If this is a linear problem (i.e. stsDMMGSetKSP() was used) then this is the 
     master linear solver. If this is a nonlinear problem (i.e. stsDMMGSetSNES() was used) this
     returns the KSP (linear solver) that is associated with the SNES (nonlinear solver)

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ(), KSPGetSNES()

M*/
#define stsDMMGGetKSP(ctx)            (ctx)[(ctx)[0]->nlevels-1]->ksp

/*MC
   stsDMMGGetSNES - Gets the SNES object (nonlinear solver) for the finest level

   Synopsis:
   SNES stsDMMGGetSNES(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Notes: If this is a linear problem (i.e. stsDMMGSetKSP() was used) then this returns PETSC_NULL

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ(), KSPGetKSP()

M*/
#define stsDMMGGetSNES(ctx)           (ctx)[(ctx)[0]->nlevels-1]->snes

/*MC
   stsDMMGGetDA - Gets the DA object on the finest level

   Synopsis:
   DA stsDMMGGetDA(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Notes: Use only if the stsDMMG was created with a DA, not a VecPack

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ(), KSPGetKSP(), stsDMMGGetVecPack()

M*/
#define stsDMMGGetDA(ctx)             (DA)((ctx)[(ctx)[0]->nlevels-1]->dm)

/*MC
   stsDMMGGetVecPack - Gets the VecPack object on the finest level

   Synopsis:
   VecPack stsDMMGGetVecPack(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

   Notes: Use only if the stsDMMG was created with a DA, not a VecPack

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetJ(), KSPGetKSP(), stsDMMGGetDA()

M*/
#define stsDMMGGetVecPack(ctx)        (VecPack)((ctx)[(ctx)[0]->nlevels-1]->dm)

/*MC
   stsDMMGGetUser - Returns the user context for a particular level

   Synopsis:
   void* stsDMMGGetUser(stsDMMG *dmmg,PetscInt level)

   Not Collective

   Input Parameters:
+   dmmg - stsDMMG solve context
-   level - the number of the level you want the context for

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser()

M*/
#define stsDMMGGetUser(ctx,level)     ((ctx)[level]->user)

/*MC
   stsDMMGSetUser - Sets the user context for a particular level

   Synopsis:
   PetscErrorCode stsDMMGSetUser(stsDMMG *dmmg,PetscInt level,void *ctx)

   Not Collective

   Input Parameters:
+   dmmg - stsDMMG solve context
.   level - the number of the level you want the context for
-   ctx - the context

   Level: intermediate

   Note: if the context is the same for each level just pass it in with 
         stsDMMGCreate() and don't call this macro

.seealso: stsDMMGCreate(), stsDMMGGetUser()

M*/
#define stsDMMGSetUser(ctx,level,usr) ((ctx)[level]->user = usr,0)

/*MC
   stsDMMGGetLevels - Gets the number of levels in a stsDMMG object

   Synopsis:
   PetscInt stsDMMGGetLevels(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGGetUser()

M*/
#define stsDMMGGetLevels(ctx)         (ctx)[0]->nlevels

/*MC
   stsDMMGGetstsDMMG - Returns the stsDMMG struct for the finest level

   Synopsis:
   stsDMMG stsDMMGGetstsDMMG(stsDMMG *dmmg)

   Not Collective

   Input Parameter:
.   dmmg - stsDMMG solve context

   Level: intermediate

.seealso: stsDMMGCreate(), stsDMMGSetUser(), stsDMMGGetB()

M*/
#define stsDMMGGetstsDMMG(ctx)              (ctx)[(ctx)[0]->nlevels-1]

PETSC_EXTERN_CXX_END
#endif
