TARGET = markTags
version += 0.1
TEMPLATE =   app
CONFIG += release

FORMS += dcmDlg.ui \
          saveVolumeDlg.ui 
SOURCES += main.cpp \
           xform.cpp \
           arthurstyle.cpp \
           arthurwidgets.cpp \
           hoverpoints.cpp \
           Volume.cpp \
           dbh.cpp \ 
           Point.cpp \
           dcmDialog.cpp \
           saveVolumeDialog.cpp

HEADERS += xform.h \
           arthurstyle.h \
           arthurwidgets.h \
           hoverpoints.h \
           Volume.h \
           dbh.h \
           Point.h \
           dcmDialog.h \
           saveVolumeDialog.h \
           dcmDataSeries.h

RESOURCES += shared.qrc

config = 
unix {
  DEFINES =     HAVE_CONFIG_H
  LIBPATH =     /usr/local/lib \
    /usr/local/dicom/lib
  LIBS = -ldcmdata \
         -lofstd
  INCLUDEPATH = /usr/local/dicom/include
}
win32 {
  DEFINES = _WIN32
  LIBPATH =  $(DCMTK_DIR)/lib
  LIBS =    -ldcmdata \
          -lofstd \
          -lnetapi32 \
        -lwsock32
  INCLUDEPATH =  $(DCMTK_DIR)/include
}

QT +=   core \
  xml \
  gui \
  opengl

