#include "tooldock.h"
#include "viewer3d.h"


toolDock::toolDock( const QString &name, QWidget * parent, Qt::WFlags f) : QWidget(parent, f) //QDockWidget(parent, f)
{
 	setObjectName(name);
 	// setWindowTitle(name);

	setupUi(this);
}

void toolDock::setVolDims(QString dims) {
	dimensions->setText(dims);
}

void toolDock::setVolSpacing(QString sp) {
	pixSpacing->setText(sp);
}

void toolDock::setPatientName(QString pn) {
	patientName->setText(pn);
}

void toolDock::setMaxDepth(int i) {
  octLowSpinBox->setMaximum(i);
  octHighSpinBox->setMaximum(i);
  
  horizontalSlider_2->setMaximum(i);
  horizontalSlider_3->setMaximum(i);  
  
  octHighSpinBox->setValue(i);
}

void toolDock::showTab(int i) {
  methodsTabWidget->setCurrentIndex(i);
}

//
void toolDock::on_fibersCheckBox_toggled(bool)
{
}
//
void toolDock::on_jacCheckBox_toggled(bool )
{
}
//
void toolDock::on_defCheckBox_toggled(bool )
{
}
//
void toolDock::on_clipCheckBox_toggled(bool )
{
	m_viewer->toggleClipPlane();
}
//
void toolDock::on_bboxCheckBox_toggled(bool)
{
	m_viewer->toggleBoundingBox();
}
//
void toolDock::on_isoCheckBox_toggled(bool )
{
}
//
void toolDock::on_volCheckBox_toggled(bool)
{
	   m_viewer->toggleVolume();
}
//
void toolDock::on_octCheckBox_toggled(bool)
{
	   m_viewer->toggleOctree();
}
//
void toolDock::on_textBrowser_anchorClicked(QUrl )
{
}
//
void toolDock::on_configServersButton_clicked()
{
}
//
void toolDock::on_submitJobButton_clicked()
{
}
//
void toolDock::on_viewCompletedJobsButton_clicked()
{
}
//
void toolDock::on_renOctListWidget_itemDoubleClicked(QListWidgetItem* item)
{
}
//
void toolDock::on_renOctListWidget_itemSelectionChanged()
{
}
//
void toolDock::on_octLowSpinBox_valueChanged(int )
{
  // update the viewer
  m_viewer->setOctThreshold( octLowSpinBox->value(), octHighSpinBox->value() );
}
//
void toolDock::on_octHighSpinBox_valueChanged(int )
{
  // update the viewer
  m_viewer->setOctThreshold( octLowSpinBox->value(), octHighSpinBox->value() );
}
//
void toolDock::on_seqRewButton_clicked()
{
	m_viewer->controlSequencer(viewer3d::REW);
}
//
void toolDock::on_seqPlayButton_clicked()
{
	m_viewer->controlSequencer(viewer3d::PLAY);
}
//
void toolDock::on_seqPauseButton_clicked()
{
	m_viewer->controlSequencer(viewer3d::PAUSE);
}
//
void toolDock::on_seqFwdButton_clicked()
{
	m_viewer->controlSequencer(viewer3d::FWD);
}
//
void toolDock::on_frametRateSpinBox_valueChanged(double )
{
}
//
void toolDock::on_isoApplyAllButton_clicked()
{
    m_viewer->updateIsoAll(isoSpinBox->value());
}
//
void toolDock::on_isoOpacitySpinBox_valueChanged(double )
{
}
//
void toolDock::on_isoSpinBox_valueChanged(int iso)
{
	m_viewer->changeIsoLevel(iso);
}
//
void toolDock::on_addIsoButton_clicked()
{
}
//
void toolDock::on_loadPDButton_clicked()
{
	QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "", "MetaIO Files (*.mhd)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadPD( fn );
        fiberTrackButton->setEnabled(true);
	}
}
//
void toolDock::on_loadFAButton_clicked()
{
	QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "", "MetaIO Files (*.mhd)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadFA( fn );
	}
}
//
void toolDock::on_processDTButton_clicked()
{
}
//
void toolDock::on_loadDTButton_clicked()
{
  QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "", "Point files (*.pts)" );
  if ( !fn.isEmpty() ) {
    m_viewer->loadPoints( fn );
  }
}
//
void toolDock::on_fiberTrackButton_clicked()
{
	m_viewer->fiberTrack(fiberSamplingSpinBox->value(), fiberColorComboBox->currentIndex());
}
//
void toolDock::on_estimateFibersButton_clicked()
{
}
//
void toolDock::on_listWidget_itemDoubleClicked(QListWidgetItem* item)
{
}
//
void toolDock::on_presetSaveButton_clicked()
{
}
//
void toolDock::on_deletePresetButton_clicked()
{
}
//
void toolDock::on_tfeResetButton_clicked()
{
	gradienteditor->reset();
}
//
void toolDock::on_gradienteditor_colormapChanged()
{
  // get colormap pointer from the viewer
  unsigned char * cmap = m_viewer->getColormap(); 
  // update it via the tfEditor.
  gradienteditor->updateColormap(cmap, 256); // size is hardcoded for now ... will allways be RGBA.

  m_viewer->colormapChanged();
}
//

void toolDock::on_defVisSparsitySlider_valueChanged(int) {
	m_viewer->redrawDeformationField();
}

void toolDock::on_defVisSliceSlider_valueChanged(int) {
	m_viewer->redrawDeformationField();
}

void toolDock::on_defVisPlaneCombo_currentIndexChanged(int) {
	m_viewer->redrawDeformationField();
}

void toolDock::on_defVisModeCombo_currentIndexChanged(int){
	m_viewer->redrawDeformationField();
}

void toolDock::on_defVisColorCombo_currentIndexChanged(int) {
	m_viewer->redrawDeformationField();
}

void toolDock::on_defFieldCombo_currentIndexChanged(int) {
	m_viewer->redrawDeformationField();
}

void toolDock::on_defVisFieldNext_clicked() {
	int val = defFieldCombo->currentIndex() +1;
	if (val >= defFieldCombo->count() )
		val =0;
	defFieldCombo->setCurrentIndex( val );
}

void toolDock::on_defVisFieldPrev_clicked() {
	int val = defFieldCombo->currentIndex() -1;
	if (val < 0 )
		val = defFieldCombo->count() -1;
	defFieldCombo->setCurrentIndex( val );
}