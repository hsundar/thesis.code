#ifndef PATIENTBROWSER_H
#define PATIENTBROWSER_H
//
#include "ui_patientBrowser.h"
#include "dcmDataSeries.h"
#include <QDir>
//
class patientBrowser : public QDialog, public Ui::patientBrowser
{
Q_OBJECT
public:
	patientBrowser( QWidget * parent = 0, Qt::WFlags f = 0 );
	bool initDataDirectory(QDir dir);
private slots:
	void on_deleteButton_clicked();
	void on_importButton_clicked();
	void on_cdImportButton_clicked();
	void on_editPatientButton_clicked();
	void on_viewButton_clicked();
	void on_reloadButton_clicked();
	void on_databaseButton_clicked();
	void on_exportButton_clicked();
	void on_anonymizeButton_clicked();
	void on_hdImportButton_clicked();
        
        void resizeCols();
        void updatePreview();

protected:
	bool rescanDataDirectory();
	
  // variables ......
	QDir						directory;
	QList<dcmDataSeries> dseries;
private:
	void addFilesFromDir( const QString& directory, QStringList& fileList);
	void scanDir(const QString& directory, QStringList& fileList); 
};
#endif
