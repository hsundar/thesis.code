#ifndef __SAVE_VOLUME_DIALOG_H_
#define __SAVE_VOLUME_DIALOG_H_

#include "ui_saveVolumeDlg.h"

class QDirModel;

class saveVolumeDialog : public QDialog {
  Q_OBJECT

  protected slots:
    void createDirectory();
    void remove();
  
  public:
        saveVolumeDialog(QWidget *parent = 0);
        QString getSaveFormat();
        QString getSavePath();
        QString getFilePrefix();

  private:
        QDirModel *m_model;
        Ui::saveVolumeDlg ui;
};

#endif
