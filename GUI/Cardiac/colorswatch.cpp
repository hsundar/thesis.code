#include "colorswatch.h"
#include "tfe.h"
#include "viewer3d.h"
#include "mainwindow.h"

#include <QtGui>
#include <fstream>
/*
#include <QAction>
#include <QtEvents>
#include <QFrame>
#include <QMainWindow>
#include <QMenu>
#include <QPainter>
#include <QImage>
#include <QColor>
#include <QPainterPath>
#include <QtDebug>
*/

  ColorSwatch::ColorSwatch(const QString &dockName, QWidget *parent, Qt::WFlags flags)
: QDockWidget(parent, flags)
{
  setObjectName(dockName);
  setWindowTitle(objectName());

  MainWindow *mainWin = (MainWindow *)parent;
  m_viewer = mainWin->getViewer();

  QFrame * frame = new QFrame(this);
  frame->setMinimumWidth(200);
  QVBoxLayout *layout = new QVBoxLayout( frame );
  layout->setSpacing(1);
  layout->setMargin(1);

  // INFO
  QGroupBox *infoGroup = new QGroupBox( "Information", frame);
  infoGroup->setAttribute(Qt::WA_ContentsPropagated);

  m_lblPatientName = new QLabel("<b>Patient Name:</b>", infoGroup);
  m_lblVolDims = new QLabel("<b>Dimensions:</b>", infoGroup);
  m_lblVolSpacing = new QLabel("<b>Spacing:</b>", infoGroup);

  QVBoxLayout *infoLayout = new QVBoxLayout( infoGroup );
  infoLayout->addWidget(m_lblPatientName);
  infoLayout->addWidget(m_lblVolDims);
  infoLayout->addWidget(m_lblVolSpacing);

  layout->addWidget(infoGroup);

  // DISPLAY
  QGroupBox *displayGroup = new QGroupBox("Display Options", frame);
  displayGroup->setAttribute(Qt::WA_ContentsPropagated);

  m_chkBBox = new QCheckBox( "Bounding Box", displayGroup );
  m_chkBBox->resize( m_chkBBox->sizeHint() );
  connect( m_chkBBox, SIGNAL (clicked() ), m_viewer, SLOT( toggleBoundingBox() ) );
  m_chkBBox->setChecked (TRUE);

  m_chkClip = new QCheckBox( "Clip Plane", displayGroup );
  m_chkClip->resize( m_chkClip->sizeHint() );
  connect( m_chkClip, SIGNAL (clicked() ), m_viewer, SLOT( toggleClipPlane() ) );

  m_chkVolume = new QCheckBox( "Volume", displayGroup );
  m_chkVolume->resize(m_chkVolume->sizeHint() );
  connect( m_chkVolume, SIGNAL (clicked()), m_viewer, SLOT( toggleVolume() ) );
  m_chkVolume->setChecked(TRUE);

  m_chkOctree = new QCheckBox( "Octree", displayGroup );
  m_chkOctree->resize( m_chkOctree->sizeHint() );
  connect( m_chkOctree, SIGNAL( clicked() ), m_viewer, SLOT( toggleOctree() ) );
  m_chkOctree->setChecked(TRUE);

  m_chkIso = new QCheckBox( "Isosurface", displayGroup );
  m_chkIso->resize( m_chkIso->sizeHint() );

  m_chkDefFld = new QCheckBox( "Deformation", displayGroup );
  m_chkDefFld->resize( m_chkDefFld->sizeHint() );

  m_chkFibers = new QCheckBox( "Fibers", displayGroup );
  m_chkFibers->resize( m_chkFibers->sizeHint() );

  m_chkJac = new QCheckBox( "|Jacobian|", displayGroup );
  m_chkJac->resize( m_chkJac->sizeHint() );

  QGridLayout *displayGroupLayout = new QGridLayout( displayGroup );
  displayGroupLayout->addWidget(m_chkBBox, 0, 0);
  displayGroupLayout->addWidget(m_chkClip, 1, 0);
  displayGroupLayout->addWidget(m_chkVolume, 2, 0);
  displayGroupLayout->addWidget(m_chkOctree, 3, 0);
  displayGroupLayout->addWidget(m_chkIso, 0, 1);
  displayGroupLayout->addWidget(m_chkDefFld, 1, 1);
  displayGroupLayout->addWidget(m_chkJac, 2, 1);
  displayGroupLayout->addWidget(m_chkFibers, 3, 1);

  layout->addWidget(displayGroup);

  // TABS
  m_tabWidget = new QTabWidget (frame);
  // m_tabWidget = new QToolBox (frame);
  // m_tabWidget->setStyle(QStyleFactory::create("Plastique"));
  // m_tabWidget->setStyle(QStyleFactory::create("WindowsXP"));

  QWidget *volTab = new QWidget (m_tabWidget);
  QWidget *octTab = new QWidget (m_tabWidget);
  QWidget *dtiTab = new QWidget (m_tabWidget);

  // Vol Tab
  QVBoxLayout *vtl = new QVBoxLayout(volTab);

  // isosurface ...
  QGroupBox *isoSurfGp = new QGroupBox("Isosurface", volTab);
  QGridLayout *isolay = new QGridLayout(isoSurfGp);

  m_sldrIso = new QSlider( Qt::Horizontal, isoSurfGp);
  m_sldrIso->setTickInterval(5);
  m_sldrIso->setMinimum(0); m_sldrIso->setMaximum(256);
  m_sldrIso->setValue(100);

  QRadioButton *m_isoApplyAllBut = new QRadioButton("Apply to all frames", isoSurfGp);
  m_isoApplyAllBut->setChecked(false);

  connect( m_isoApplyAllBut, SIGNAL(toggled(bool)), this, SLOT(isoApplyAllToggled(bool)) );


  connect( m_sldrIso, SIGNAL( valueChanged(int) ), this, SLOT (sldrIsoChanged(int)) );
  
  m_lblIso = new QLabel("100", isoSurfGp);

  isolay->addWidget(m_sldrIso, 0, 0);
  isolay->addWidget(m_lblIso, 0, 1);
  isolay->addWidget(m_isoApplyAllBut, 1, 0);

  QGroupBox *volSeqGp = new QGroupBox("Sequence", volTab);
  volSeqGp->setAttribute(Qt::WA_ContentsPropagated);
  QHBoxLayout *seqGpLayout = new QHBoxLayout(volSeqGp);

  // add 4 buttons to the group ...
  QPushButton *rewBut =  new QPushButton(QIcon(QPixmap("images/player_rew.png")), "", volSeqGp);
  QPushButton *playBut =  new QPushButton(QIcon(QPixmap("images/player_play.png")), "", volSeqGp);
  QPushButton *pauseBut =  new QPushButton(QIcon(QPixmap("images/player_pause.png")), "", volSeqGp);
  QPushButton *fwdBut =  new QPushButton(QIcon(QPixmap("images/player_fwd.png")), "", volSeqGp);
  
  // Button Options ...
  rewBut->setFlat(true); // rewBut->setIconSize(QSize(24,24));
  playBut->setFlat(true); // playBut->setIconSize(QSize(24,24));  
  pauseBut->setFlat(true); // pauseBut->setIconSize(QSize(24,24));
  fwdBut->setFlat(true); // fwdBut->setIconSize(QSize(24,24));

  // connections ...
  connect( playBut, SIGNAL( clicked() ), m_viewer, SLOT(playSeq()) );
  connect( pauseBut, SIGNAL( clicked() ), m_viewer, SLOT(pauseSeq()) );
  connect( rewBut, SIGNAL( clicked() ), m_viewer, SLOT( rewSeq() ));
  connect( fwdBut, SIGNAL( clicked() ), m_viewer, SLOT( fwdSeq()));

  vtl->addWidget(isoSurfGp);
  vtl->addWidget(volSeqGp);

  // layout ...
  seqGpLayout->addWidget(rewBut);
  seqGpLayout->addWidget(playBut);
  seqGpLayout->addWidget(pauseBut);
  seqGpLayout->addWidget(fwdBut);
  // Oct Tab

  QGroupBox *octThrGp = new QGroupBox("Octree Threshold", octTab);
  octThrGp->setAttribute(Qt::WA_ContentsPropagated);
  QVBoxLayout *thrGpLayout = new QVBoxLayout(octTab);

  QGridLayout *otl = new QGridLayout(octThrGp);
  m_sldrOctHigh = new QSlider( Qt::Horizontal, octThrGp);
  m_sldrOctHigh->setTickInterval(1);
  m_sldrOctHigh->setMinimum(0); m_sldrOctHigh->setMaximum(8);
  m_sldrOctHigh->setValue(8);
  m_sldrOctLow  = new QSlider( Qt::Horizontal, octThrGp);
  m_sldrOctLow->setTickInterval(1);
  m_sldrOctLow->setMinimum(0); m_sldrOctLow->setMaximum(8);

  connect( m_sldrOctLow, SIGNAL( valueChanged(int) ), this, SLOT (sldrLowChanged(int)) );
  connect( m_sldrOctHigh, SIGNAL( valueChanged(int) ), this, SLOT (sldrHighChanged(int)) );

  m_lblOctHigh = new QLabel("8", octThrGp);
  m_lblOctLow = new QLabel("0", octThrGp);

  otl->addWidget(new QLabel("Low", octThrGp), 0, 0);
  otl->addWidget(new QLabel("High", octThrGp), 1, 0);
  otl->addWidget(m_sldrOctLow, 0, 1);
  otl->addWidget(m_lblOctLow, 0, 2);
  otl->addWidget(m_sldrOctHigh, 1, 1);
  otl->addWidget(m_lblOctHigh, 1, 2);
  otl->setColumnStretch(1, 4);

  thrGpLayout->addWidget(octThrGp);
  // DTI Tab

  QGridLayout *dtl = new QGridLayout(dtiTab);
  QPushButton *pdBut = new QPushButton("Load PD", dtiTab);
  QPushButton *faBut = new QPushButton("Load FA", dtiTab);
  QPushButton *dtBut = new QPushButton("Load Tensors", dtiTab);
  QPushButton *trackBut = new QPushButton("Fiber Track", dtiTab);

  connect( faBut, SIGNAL(clicked()), this, SLOT(loadFA()) );
  connect( pdBut, SIGNAL(clicked()), this, SLOT(loadPD()) );
  connect( trackBut, SIGNAL(clicked()), this, SLOT(fiberTrack()) );

  m_sldrFiberDensity = new QSlider( Qt::Horizontal, dtiTab );
  m_sldrFiberDensity->setMinimum(1); m_sldrFiberDensity->setMaximum(10);
  m_sldrFiberDensity->setTickInterval(1);

  QSpinBox *fiberDensitySpin = new QSpinBox(dtiTab);
  fiberDensitySpin->setRange(1,10);

  connect(fiberDensitySpin, SIGNAL(valueChanged(int)), m_sldrFiberDensity, SLOT(setValue(int))); 
  connect(m_sldrFiberDensity, SIGNAL(valueChanged(int)), fiberDensitySpin, SLOT(setValue(int))); 
  
  m_sldrFiberDensity->setValue(4);

  m_FiberColorCombo = new QComboBox(dtiTab);
  m_FiberColorCombo->insertItem(0, "Gray");
  m_FiberColorCombo->insertItem(1, "Normal");
  m_FiberColorCombo->insertItem(2, "Chiral");
  m_FiberColorCombo->setCurrentIndex(1);


  dtl->addWidget(pdBut, 0, 0, 1, 2);
  dtl->addWidget(faBut, 0, 2, 1, 2);
  dtl->addWidget(dtBut, 1, 0, 1, 2);
  dtl->addWidget(m_FiberColorCombo, 1, 2, 1, 2);
  dtl->addWidget(m_sldrFiberDensity, 2, 0, 1, 3);
  dtl->addWidget(fiberDensitySpin, 2, 3);
  dtl->addWidget(trackBut, 3, 0, 1, 3);

  // ----
  //m_tabWidget->addItem(volTab, "&Volume");
  m_tabWidget->addTab(volTab, "&Volume");
  //m_tabWidget->addItem(octTab, "&Octree");
  //m_tabWidget->addItem(dtiTab, "&DTI" );
  m_tabWidget->addTab(octTab, "&Octree");
  m_tabWidget->addTab(dtiTab, "&DTI" );

  layout->addWidget(m_tabWidget);

  // TRANSFER FUNCTION
  QGroupBox *tfGroup = new QGroupBox("Transfer Function", frame);
  tfGroup->setAttribute(Qt::WA_ContentsPropagated);

  m_tfe = new tfEditor ( tfGroup );
  // swatch->setFrameStyle(QFrame::Box | QFrame::Sunken);
  // m_tfe->setMinimumSize(200, 200);

  connect( m_tfe, SIGNAL ( colormapChanged() ), this, SLOT (tfChanged()) ); 

  QVBoxLayout *tfGroupLayout = new QVBoxLayout( tfGroup );
  tfGroupLayout->addWidget( m_tfe );

  QPushButton* tfReset = new QPushButton("Reset", tfGroup);
  tfGroupLayout->addWidget( tfReset );

  QComboBox *tfPresets = new QComboBox( tfGroup );
  tfPresets->addItem("Cardiac MR");
  tfPresets->addItem("Cardiac CT");
  tfGroupLayout->addWidget( tfPresets );

  connect( tfReset, SIGNAL( clicked() ), m_tfe, SLOT( reset() ) );

  emit tfReset->animateClick();

  layout->addWidget(tfGroup);

  setWidget(frame);
}

void ColorSwatch::setVolDims(QString dims) {
  m_lblVolDims->setText(dims);
}

void ColorSwatch::setVolSpacing(QString sp) {
  m_lblVolSpacing->setText(sp);
}

void ColorSwatch::setPatientName(QString pn) {
  m_lblPatientName->setText(pn);
}

void ColorSwatch::sldrLowChanged(int l) {
  // update the label ...
  m_lblOctLow->setText( QString::number(l) );
  // and the viewer
  m_viewer->setOctThreshold( m_sldrOctLow->value(), m_sldrOctHigh->value() );
}

void ColorSwatch::sldrHighChanged(int h) {
  // update the label ...
  m_lblOctHigh->setText( QString::number( h ) );
  // and the viewer
  m_viewer->setOctThreshold( m_sldrOctLow->value(), m_sldrOctHigh->value() );
}

void ColorSwatch::sldrIsoChanged(int iso) {
  // update the label ...
  m_lblIso->setText( QString::number( iso ) );
  // and the viewer ... who has the isosurface ...
  m_viewer->changeIsoLevel(iso);
}

void ColorSwatch::tfChanged() {
  // get colormap pointer from the viewer
  unsigned char * cmap = m_viewer->getColormap(); 
  // update it via the tfEditor.
  m_tfe->updateColormap(cmap, 256); // size is hardcoded for now ... will allways be RGBA.

  m_viewer->colormapChanged();
}

void ColorSwatch::showTab(int i) {
  // QMessageBox::information( this, "Cardiac Viewer", "About to change tab" );
  m_tabWidget->setCurrentIndex(i);
}

void ColorSwatch::setMaxDepth(int i) {
  m_sldrOctHigh->setMaximum(i);
  m_sldrOctLow->setMaximum(i);
  
  m_lblOctLow->setText("0");
  m_lblOctHigh->setText( QString::number(i) ); 
}

void ColorSwatch::isoApplyAllToggled(bool flag) {
  if (flag)
    m_viewer->updateIsoAll(m_sldrIso->value());
}

void ColorSwatch::loadFA() {
	QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "", "MetaIO Files (*.mhd)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadFA( fn );
	}
}

void ColorSwatch::loadPD() {
	QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "", "MetaIO Files (*.mhd)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadPD( fn );
	}
}

void ColorSwatch::fiberTrack() {
  // change to pass 2 additional args ...
  // density and color ...

  m_viewer->fiberTrack(m_sldrFiberDensity->value(), m_FiberColorCombo->currentIndex());
}
