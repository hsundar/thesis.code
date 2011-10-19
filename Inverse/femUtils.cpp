#include "femUtils.h"

int elementToNode( DA da, Vec elementVec, Vec nodeVec) {
#ifdef __DEBUG__  
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif  
  int ierr;
  int x,y,z,m,n,p;
  int mx,my,mz, xne, yne, zne;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int dof=1;
  PetscScalar ***elem;
  PetscScalar ***node;
  CHKERRQ( DAGetCorners(da, &x, &y, &z, &m, &n, &p) ); 
  CHKERRQ( DAGetInfo(da,0, &mx, &my, &mz, 0,0,0,&dof,0,0,0) ); 
  if (x+m == mx) {
    xne=m-1;
  } else {
    xne=m;
  }
  if (y+n == my) {
    yne=n-1;
  } else {
    yne=n;
  }
  if (z+p == mz) {
    zne=p-1;
  } else {
    zne=p;
  }
  Vec localVec;

  CHKERRQ( DACreateLocalVector(da, &localVec) );
  CHKERRQ( VecSet( localVec, 0.0));
  ierr = DAVecGetArray(da, elementVec,&elem); CHKERRQ(ierr);
  ierr = DAVecGetArray(da, localVec,&node); CHKERRQ(ierr);
  for (int k=z; k<z+zne; k++) {
    for (int j=y; j<y+yne; j++) {
      for (int i=x; i<x+xne; i++) {
        for (int d=0; d<dof; d++) {
          node[k][j][dof*i+d]         += elem[k][j][dof*i+d]/8.0; 
          node[k][j][dof*(i+1)+d]     += elem[k][j][dof*i+d]/8.0;
          node[k][j+1][dof*i+d]       += elem[k][j][dof*i+d]/8.0; 
          node[k][j+1][dof*(i+1)+d]   += elem[k][j][dof*i+d]/8.0;
          node[k+1][j][dof*i+d]       += elem[k][j][dof*i+d]/8.0; 
          node[k+1][j][dof*(i+1)+d]   += elem[k][j][dof*i+d]/8.0;
          node[k+1][j+1][dof*i+d]     += elem[k][j][dof*i+d]/8.0; 
          node[k+1][j+1][dof*(i+1)+d] += elem[k][j][dof*i+d]/8.0;
        } // d
      }  // i
    } // j
  } // k
  
  ierr = DAVecRestoreArray(da, elementVec, &elem); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da, localVec, &node); CHKERRQ(ierr);
  ierr = DALocalToGlobalBegin(da, localVec, nodeVec); CHKERRQ(ierr);
  ierr = DALocalToGlobalEnd(da, localVec, nodeVec); CHKERRQ(ierr);
#ifdef __DEBUG__  
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif  
  return(0);
}

int elementToNode(ot::DA da, Vec elementVec, Vec nodeVec, unsigned int dof) {
#ifdef __DEBUG__  
  std::cout << RED"Entering "NRM << __func__ << std::endl;
#endif  
  PetscScalar *elem, *node;

  da.vecGetBuffer(elementVec, elem, true, false, true, dof);
  da.vecGetBuffer(nodeVec, node, false, true, false, dof);

  // loop ...
  for ( da.init<ot::DA::ALL>(), da.init<ot::DA::WRITABLE>(); da.curr() < da.end<ot::DA::ALL>(); da.next<ot::DA::ALL>()) {
    unsigned int currIndex = da.curr();
    unsigned int indices[8];
    da.getNodeIndices(indices); 
    
    unsigned char hn = da.getHangingNodeIndex(da.curr());

    for (int d=0; d<dof; d++) {
      for(int i = 0; i < 8; i++) {
        if (!(hn & (1 << i))) {          
            node[dof*indices[i]+d] += elem[dof*currIndex + d]/8;
        }
      }
    }
  }

  // write to ghost nodes
  da.WriteToGhostsBegin(node, dof);
  da.WriteToGhostsEnd(node, dof);

  da.vecRestoreBuffer(elementVec, elem, true, true, true, dof);
  da.vecRestoreBuffer(nodeVec, node, false, true, false, dof);

#ifdef __DEBUG__  
  std::cout << GRN"Leaving "NRM << __func__ << std::endl;
#endif  
  return(0);
}

int nodeToElement(DA da, Vec nodeVec, Vec elementVec) {
  return 0;
}

/*
 * Functions to read and write C-array like structures in parallel. These functions 
 * are useful when all processors want to read from the same file. 
 * 
 * Move these functions over to parUtils eventually.
 */
int readParallel(DA da, char *filename, Vec vec, MPI_Comm comm) {
  // Root will read the file and scatter it.
  return 0;
}

int writeParallel(DA da, char *filename, Vec vec, MPI_Comm comm) {
  return 0;
}

int readParallel(ot::DA da, char *filename, Vec vec, MPI_Comm comm) {
  return 0;
}

int writeParallel(ot::DA da, char *filename, Vec vec, MPI_Comm comm) {
  return 0;
}

template <typename T>
int readParallel(ot::DA da, char *filename, std::vector<T> &vec, MPI_Comm comm) {
  return 0;
}

template <typename T>
int writeParallel(ot::DA da, char *filename, std::vector<T> &vec, MPI_Comm comm) {
  return 0;
}

int concatenateVecs(std::vector<Vec> dyna, Vec &ctrl, bool createNew) {
#ifdef __DEBUG__  
  std::cout << "Entering " << __func__ << std::endl;
#endif
  // First get the local size of the Vec. All time points will have the same
  // partition and consequently the same size.
  PetscInt lSize=0;
  VecGetLocalSize( dyna[0], &lSize );
  
  // std::cout << YLW"dyna size is "NRM << dyna.size() << std::endl;

  // The new local size of the vector
  PetscInt newSize = dyna.size()*lSize;
  
  // std::cout << "Sizes are " << lSize << ", " << newSize << std::endl;

  // Now create the ctrl vector.
  if (createNew) {
    VecCreate(PETSC_COMM_WORLD, &ctrl);
    VecSetSizes(ctrl, newSize, PETSC_DECIDE);
    VecSetFromOptions(ctrl);
  }
  // std::cout << "Created Vector, Now copying data" << std::endl;
  // Now we copy the data to the new vector
  PetscScalar* ctrlArray, *dynaArray;
  VecGetArray(ctrl, &ctrlArray);
  for (unsigned int i=0; i<dyna.size(); i++) {
    VecGetArray(dyna[i], &dynaArray);
    for (int j=0; j<lSize; j++)
      ctrlArray[i*lSize + j] = dynaArray[j];
    VecRestoreArray(dyna[i], &dynaArray);
  }
  VecRestoreArray(ctrl, &ctrlArray);
  
  // VecGetLocalSize( ctrl, &newSize );
#ifdef __DEBUG__  
  std::cout << "Leaving " << __func__ << std::endl;
#endif
  return 0;
}

int splitVec(Vec ctrl, Vec* &dyna, unsigned int nt) {
  // Vec x;
  // clear the std::vector to be safe
  dyna = (Vec *)malloc(nt*sizeof(Vec));
  // Get the right sizes ...
  PetscInt lSize=0;
  VecGetLocalSize( ctrl, &lSize );
  // The new local size of the vector
  PetscInt newSize = lSize/nt;

#ifdef __DEBUG__
  if ( newSize*nt != lSize ) {
    std::cerr << __func__ << ": Trying to split into non-integer nodes." << std::endl;
    return 1;
  }
#endif
  
  PetscScalar* ctrlArray, *dynaArray;
  VecGetArray(ctrl, &ctrlArray);
  for (unsigned int i=0; i<nt; i++) {
    // Create the Vector ...
    VecCreate(PETSC_COMM_WORLD, &dyna[i]);
    VecSetSizes(dyna[i], newSize, PETSC_DECIDE);
    VecSetFromOptions(dyna[i]);
    VecGetArray(dyna[i], &dynaArray);
    for (int j=0; j<lSize; j++)
      dynaArray[j] = ctrlArray[i*lSize + j];
    VecRestoreArray(dyna[i], &dynaArray);
  }
  VecRestoreArray(ctrl, &ctrlArray);
 
  return 0;
}

int splitVec(Vec ctrl, std::vector<Vec> &dyna, unsigned int nt) {
  // std::cout << "Entering " << __func__ << std::endl;
  Vec x;
  // clear the std::vector to be safe
  dyna.clear();
  // Get the right sizes ...
  PetscInt lSize=0;
  VecGetLocalSize( ctrl, &lSize );
  // The new local size of the vector
  PetscInt newSize = lSize/nt;

  // std::cout << "Sizes are " << lSize << ", " << newSize << std::endl;

#ifdef __DEBUG__
  if ( newSize*nt != lSize ) {
    std::cerr << __func__ << ": Trying to split into non-integer nodes." << std::endl;
    return 1;
  }
#endif
  
  PetscScalar *ctrlArray, *dynaArray;
  VecGetArray(ctrl, &ctrlArray);
  for (unsigned int i=0; i<nt; i++) {
    // Create the Vector ...
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, newSize, PETSC_DECIDE);
    VecSetFromOptions(x);
    VecGetArray(x, &dynaArray);
    for (int j=0; j<newSize; j++)
      dynaArray[j] = ctrlArray[i*newSize + j];
    VecRestoreArray(x, &dynaArray);
    dyna.push_back(x);
  }
  VecRestoreArray(ctrl, &ctrlArray);
  // std::cout << "Leaving " << __func__ << std::endl;
 
  return 0;
}


