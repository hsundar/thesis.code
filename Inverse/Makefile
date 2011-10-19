SHELL = /bin/sh
CEXT = cpp
include ${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf
include ${PETSC_DIR}/bmake/common/variables
EXEC = genLVfibers genFiberActivation genCmameFibers genLV fwd_RG_fullForce fwd_RG_fiberForce fwd_Oct_fullForce fwd_Oct_fiberForce inv_RG_fullForce inv_RG_fiberForce inv_Oct_fullForce inv_Oct_fiberForce
CFLAGS = -O3 #-D_PETSC_USE_LOG_ #-D__DEBUG__ # -D_OCT_CHECK_ 
GC = g++
INCLUDE = -I./  -I$(OTK_DIR)/include/oct -I$(OTK_DIR)/include/stsmg -I$(OTK_DIR)/include/oda  -I$(OTK_DIR)/include/par  -I$(OTK_DIR)/include/shape  -I$(OTK_DIR)/include/petsc  -I$(OTK_DIR)/include/mat  -I$(OTK_DIR)/include/volume  -I$(OTK_DIR)/include/point  -I$(OTK_DIR)/include/test -I$(OTK_DIR)/include/binOps -I$(OTK_DIR)/include/random -I$(OTK_DIR)/include/indexHolder -I$(OTK_DIR)/include  ${PETSC_INCLUDE} #-I$(OTK_DIR)/MatVecODA
LIBS = -L$(OTK_DIR)/lib -lODA -lOct -lPar -lPoint -lTest -lBinOps -lPsc ${PETSC_LIB}

all : $(EXEC)
utils : genLVfibers genFiberActivation genCmameFibers genLV
fwd : fwd_RG_fullForce fwd_RG_fiberForce fwd_Oct_fullForce fwd_Oct_fiberForce
inv : inv_RG_fullForce inv_RG_fiberForce inv_Oct_fullForce inv_Oct_fiberForce


%.o: %.$(CEXT) 
	$(PCC) $(CFLAGS) -c $(INCLUDE) $< -o $@

# Utils
genLV : genLV.o
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@ 

genLVfibers : genLVfibers.o
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@ 

genCmameFibers : genCmameFibers.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@ 

genFiberActivation : genFiberActivation.o
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@ 

# FORWARD 
fwd_RG_fullForce : fwd_RG_fullForce.o timeStepper.o femUtils.o stsdamg.o
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@ 
	
fwd_RG_fiberForce : fwd_RG_fiberForce.o timeStepper.o femUtils.o stsdamg.o
	$(PCC) $(CFLAGS) $^ $(LIBS)	-o $@ 

fwd_Oct_fullForce : fwd_Oct_fullForce.o timeStepper.o femUtils.o stsdamg.o
	$(PCC) $(CFLAGS) $^ -o $@ $(LIBS)

fwd_Oct_fiberForce : fwd_Oct_fiberForce.o timeStepper.o femUtils.o stsdamg.o
	$(PCC) $(CFLAGS) $^ -o $@ $(LIBS)

# INVERSE
inv_RG_fullForce : inv_RG_fullForce.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@

inv_RG_fiberForce : inv_RG_fiberForce.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@

inv_Oct_fullForce : inv_Oct_fullForce.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@

inv_Oct_fiberForce : inv_Oct_fiberForce.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@

##~~~~~~~~~~

estimateCardiac : estimateCardiac.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@

linElasCheck : linElasCheck.o timeStepper.o femUtils.o stsdamg.o
	$(PCC) $(CFLAGS) $^ -o $@ $(LIBS)


elasInv : elasInverse.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ -o $@ $(LIBS)

pForce : parametricCardiac.o timeStepper.o inverseSolver.o stsdamg.o
	$(PCC) $(CFLAGS) $^ -o $@ $(LIBS)

pElas : pElas.o timeStepper.o femUtils.o inverseSolver.o stsdamg.o 
	$(PCC) $(CFLAGS) $^ $(LIBS) -o $@


clean :
	rm -rf src/*~ src/*.o *.o $(EXEC) *.exe 


