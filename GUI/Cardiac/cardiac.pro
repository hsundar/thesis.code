TARGET =   cardiac
VERSION +=   0.3
FORMS +=   dcmDlg.ui \
  saveVolumeDlg.ui \
  toolDock.ui \
  patientBrowser.ui \
  PreferencesDlg.ui 
SOURCES +=   main.cpp \
  tfe.cpp \
  colorswatch.cpp \
  mainwindow.cpp \
  Point.cpp \
  Volume.cpp \
  BoundingBox.cpp \
  MyExaminerViewer.cpp \
  viewer3d.cpp \
  arthurstyle.cpp \
  arthurwidgets.cpp \
  hoverpoints.cpp \
  dcmDialog.cpp \
  saveVolumeDialog.cpp \
  dbh.cpp \
  MarchingCubes.cpp \
  SFScalarField.cpp \
  SFVectorField.cpp \
  SoVectorFieldViz.cpp \
  plugindialog.cpp \
  FiberTracker.cpp \
  tooldock.cpp \
  ColorButton.cpp \
  patienttree.cpp \
  patientBrowser.cpp 
HEADERS +=   tfe.h \
  mainwindow.h \
  colorswatch.h \
  Point.h \
  Volume.h \
  MyExaminerViewer.h \
  BoundingBox.h \
  viewer3d.h \
  arthurstyle.h \
  arthurwidgets.h \
  hoverpoints.h \
  dbh.h \
  dcmDialog.h \
  saveVolumeDialog.h \
  MarchingCubes.h \
  SFScalarField.h \
  SFVectorField.h \
  SoVectorFieldViz.h \
  plugindialog.h \
  interfaces.h \
  Field.h \
  FiberTracker.h \
  tooldock.h \
  ColorButton.h \
  patienttree.h \
  patientBrowser.h 
config += 
unix {
  DEFINES +=     HAVE_CONFIG_H
  LIBPATH +=     /usr/local/lib \
    /usr/local/dicom/lib 
  LIBS +=     -lSimVoleon \
    -lCoin \
    -lSoQt \
    -ldcmdata \
    -lofstd 
  INCLUDEPATH +=     /usr/local/include/Inventor/annex \
    /usr/local/dicom/include
}
win32 {
  DEFINES +=     COIN_DLL \
    SOQT_DLL \
    SIMVOLEON_DLL \
    _WIN32
  LIBPATH +=     $(COINDIR)/lib \
    $(DCMTK_DIR)/lib
  LIBS +=     -lSimVoleon2 \
    -lCoin2 \
    -lSoQt1 \
    -ldcmdata \
    -lofstd \
    -lnetapi32 \
    -lwsock32
  INCLUDEPATH +=     $(COINDIR)/include \
    $(DCMTK_DIR)/include
}
RESOURCES +=   cardiac.qrc \
  shared.qrc
RC_FILE +=   icon.rc
target.path +=   $$[QT_INSTALL_DEMOS]/cardiac
sources.files +=   $$SOURCES \
  $$HEADERS \
  $$RESOURCES \
  *.pro \
  *.html
sources.path +=   $$[QT_INSTALL_DEMOS]/cardiac
INSTALLS +=   target \
  sources
CONFIG += release
QT +=   core \
  xml \
  gui \
  opengl
TEMPLATE =   app
