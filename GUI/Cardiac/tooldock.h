#ifndef TOOLDOCK_H
#define TOOLDOCK_H
//
#include "ui_toolDock.h"
//

class viewer3d;

class toolDock : public QWidget, public Ui::toolDock
{
Q_OBJECT
public:
	toolDock(const QString &name,  QWidget * parent = 0, Qt::WFlags f = 0 );
	void setViewer (viewer3d *view) { m_viewer = view; };

    // Functions for info to be set from the Viewer side ...
    void setVolDims(QString dims);
    void setVolSpacing(QString sp);
    void setPatientName(QString pn);

    void showTab(int i);
    void setMaxDepth(int i);
	
private slots:
	void on_fibersCheckBox_toggled(bool checked);
	void on_jacCheckBox_toggled(bool checked);
	void on_defCheckBox_toggled(bool checked);
	void on_clipCheckBox_toggled(bool checked);
	void on_bboxCheckBox_toggled(bool checked);
	void on_isoCheckBox_toggled(bool checked);
	void on_volCheckBox_toggled(bool checked);
	void on_octCheckBox_toggled(bool checked);
	void on_textBrowser_anchorClicked(QUrl );
	void on_configServersButton_clicked();
	void on_submitJobButton_clicked();
	void on_viewCompletedJobsButton_clicked();
	void on_renOctListWidget_itemDoubleClicked(QListWidgetItem* item);
	void on_renOctListWidget_itemSelectionChanged();
	void on_octLowSpinBox_valueChanged(int );
	void on_octHighSpinBox_valueChanged(int );
	void on_seqRewButton_clicked();
	void on_seqPlayButton_clicked();
	void on_seqPauseButton_clicked();
	void on_seqFwdButton_clicked();
	void on_frametRateSpinBox_valueChanged(double );
	void on_isoApplyAllButton_clicked();
	void on_isoOpacitySpinBox_valueChanged(double );
	void on_isoSpinBox_valueChanged(int );
	void on_addIsoButton_clicked();
	void on_loadPDButton_clicked();
	void on_loadFAButton_clicked();
	void on_processDTButton_clicked();
	void on_loadDTButton_clicked();
	void on_fiberTrackButton_clicked();
	void on_estimateFibersButton_clicked();
	void on_listWidget_itemDoubleClicked(QListWidgetItem* item);
	void on_presetSaveButton_clicked();
	void on_deletePresetButton_clicked();
	void on_tfeResetButton_clicked();
	void on_gradienteditor_colormapChanged();

	void on_defFieldCombo_currentIndexChanged(int);
	void on_defVisFieldNext_clicked();
	void on_defVisFieldPrev_clicked();
	void on_defVisSparsitySlider_valueChanged(int);
	void on_defVisSliceSlider_valueChanged(int);
	void on_defVisPlaneCombo_currentIndexChanged(int);
	void on_defVisModeCombo_currentIndexChanged(int);
	void on_defVisColorCombo_currentIndexChanged(int);
private:
    viewer3d	*m_viewer;
};
#endif
