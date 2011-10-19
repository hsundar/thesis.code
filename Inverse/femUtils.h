#ifndef __FEM_UTILS_H_
#define __FEM_UTILS_H_

#include "petscda.h"
#include "oct.h"
#include "oda.h"

/* 
 * Functions for converting a field defined on the nodes to one defined 
 * on the elements, and vice versa.
 * 
 * Interfaces for both regular grid and octree based DA structures are 
 * provided.
 *  
 */

int elementToNode(DA da, Vec elementVec, Vec nodeVec);
int nodeToElement(DA da, Vec nodeVec, Vec elementVec);

int elementToNode(ot::DA da, Vec elementVec, Vec nodeVec, unsigned int dof);
int nodeToElement(ot::DA da, Vec nodeVec, Vec elementVec, unsigned int dof);

template <typename T>
int elementToNode(ot::DA da, std::vector<T> &elementVec, std::vector<T> &nodeVec, unsigned int dof);

template <typename T>
int nodeToElement(DA da, std::vector<T> &nodeVec, std::vector<T> &elementVec, unsigned int dof);

/* 
 * Functions to convert between a field defined on a regular grid and one on a given octree. All
 * functions require that a octree-based DA be specified. 
 */

int octree2rg(ot::DA da, Vec octVec, Vec rgVec, unsigned int dof);
int rg2octree(ot::DA da, Vec rgVec, Vec octVec, unsigned int dof);

template <typename T>
int octree2rg(ot::DA da, std::vector<T> &octVec, std::vector<T> &rgVec, unsigned int dof);
template <typename T>
int rg2octree(ot::DA da, std::vector<T> &rgVec, std::vector<T> &octVec, unsigned int dof);

/*
 * Functions to read and write C-array like structures in parallel. These functions 
 * are useful when all processors want to read from the same file. 
 * 
 * Move these functions over to parUtils eventually.
 */

int readParallel(DA da, char *filename, Vec vec, MPI_Comm comm);
int writeParallel(DA da, char *filename, Vec vec, MPI_Comm comm);

int readParallel(ot::DA da, char *filename, Vec vec, MPI_Comm comm);
int writeParallel(ot::DA da, char *filename, Vec vec, MPI_Comm comm);

template <typename T>
int readParallel(ot::DA da, char *filename, std::vector<T> &vec, MPI_Comm comm);
template <typename T>
int writeParallel(ot::DA da, char *filename, std::vector<T> &vec, MPI_Comm comm);


/***********************************************
 * Functions to convert between dynamic arrays * 
 * as needed by the timesteppers and a single  *
 * PETSc Vec, which is required by the inverse *
 * solver.                                     *
 ***********************************************/

int concatenateVecs(std::vector<Vec> dyna, Vec &ctrl, bool createNew=true);
int splitVec(Vec ctrl, std::vector<Vec> &dyna, unsigned int nt);
int splitVec(Vec ctrl, Vec* &dyna, unsigned int nt);

#endif
