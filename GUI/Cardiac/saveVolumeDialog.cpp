#include "saveVolumeDialog.h"

#include <QtGui>

saveVolumeDialog::saveVolumeDialog(QWidget *parent)
: QDialog(parent)
{
  ui.setupUi(this);
  
  m_model = new QDirModel;
  m_model->setFilter( QDir::Dirs | QDir::Drives | QDir::Writable );
  m_model->setReadOnly(false);
  m_model->setSorting(QDir::DirsFirst | QDir::IgnoreCase | QDir::Name);

  ui.m_pDirView->setModel(m_model);
  ui.m_pDirView->setRootIndex(m_model->index(QDir::rootPath()));

 // ui.m_pDirView->setColumnHidden(1, true); 
  ui.m_pDirView->setColumnHidden(2, true);
  ui.m_pDirView->setColumnHidden(3, true);
 //  ui.m_pDirView->header()->hide(); 

  ui.m_pDirView->header()->setStretchLastSection(true);
  ui.m_pDirView->header()->setSortIndicator(0, Qt::AscendingOrder);
  ui.m_pDirView->header()->setSortIndicatorShown(true);
  ui.m_pDirView->header()->setClickable(true);

  QModelIndex index = m_model->index(QDir::currentPath());
  // QModelIndex index = m_model->index("/home/hsundar/data");
  ui.m_pDirView->expand(index);
  ui.m_pDirView->scrollTo(index);
  ui.m_pDirView->resizeColumnToContents(0);
  ui.m_pDirView->setCurrentIndex(index);
  
  ui.m_pDirView->show();

  connect(ui.m_pCreateDirBut, SIGNAL(clicked()), this, SLOT(createDirectory()));
  connect(ui.m_pRemoveBut, SIGNAL(clicked()), this, SLOT(remove()));
}


QString saveVolumeDialog::getSaveFormat() {
  return ui.m_pFormat->currentText();
}

QString saveVolumeDialog::getSavePath() {
  return m_model->filePath(ui.m_pDirView->currentIndex());
}

QString saveVolumeDialog::getFilePrefix() {
  return ui.m_pPrefix->text();
}

void saveVolumeDialog::createDirectory()
{
  QModelIndex index =  ui.m_pDirView->currentIndex();
  if (!index.isValid())
    return;
  QString dirName = QInputDialog::getText(this,
      tr("Create Directory"),
      tr("Directory name"));
  if (!dirName.isEmpty()) {
    if (!m_model->mkdir(index, dirName).isValid())
      QMessageBox::information(this, tr("Create Directory"),
          tr("Failed to create the directory"));
  }
}

void saveVolumeDialog::remove()
{
  QModelIndex index = ui.m_pDirView->currentIndex();
  if (!index.isValid())
    return;
  bool ok;
  if (m_model->fileInfo(index).isDir()) {
    ok = m_model->rmdir(index);
  } else {
    ok = m_model->remove(index);
  }
  if (!ok)
    QMessageBox::information(this, tr("Remove"),
        tr("Failed to remove %1").arg(m_model->fileName(index)));
}
