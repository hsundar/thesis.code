#ifndef __DCM_DIALOG_H_
#define __DCM_DIALOG_H_

#include <QtGui>

#include "ui_dcmDlg.h"

#include "dcmDataSeries.h"

#define VOL_FMT_MHD 1982
#define VOL_FMT_ANALYZE 1983


typedef struct _volInfo {
  int x,y,z;
  double sx, sy, sz;
  int bs, ba;
} volInfo;

class dcmDialog : public QDialog {
  Q_OBJECT

  protected slots:
        void changeDir();
        void scanDir();
		bool scanDirRecursive();
		void save3D();
        void save4D();
        void resizeCols();
        void updatePreview();

  public:
        dcmDialog(QWidget *parent = 0);

  protected:
        // functions ...
        void saveHeader(int format, QString fname, volInfo _info);
      
        // variables ...

        QList<dcmDataSeries> dseries;
        QDir dcmDir;
  private:
        Ui::dcmDlg ui;

		void addFilesFromDir( const QString& dcmDir, QStringList& fileList );
		void scanDir(const QString& dcmDir, QStringList& fileList);
        
};

#endif
