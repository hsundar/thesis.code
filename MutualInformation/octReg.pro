######################################################################
# Automatically generated by qmake (2.00a) Sat Feb 3 14:09:09 2007
######################################################################

TEMPLATE = app
TARGET += 
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += oct_sim.h \
           Point.h \
           TransMatrix.h \
           Volume.h \
           dbh.h \
           octThread.h \
           OptimizerPowellBrent.h \
           Optimizer.h \
           OptimizerGradientDescent.h \
           OptimizerOctSim.h \
           colors.h
SOURCES += octReg.cpp \
           oct_sim.cpp \
           Point.cpp \
           TransMatrix.cpp \
           Volume.cpp \
           dbh.cpp \
           octThread.cpp \
           OptimizerPowellBrent.cpp \
           OptimizerGradientDescent.cpp \
           Optimizer.cpp

CONFIG += release
CONFIG += thread

