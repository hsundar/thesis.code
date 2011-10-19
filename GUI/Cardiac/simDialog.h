#ifndef SIMDIALOG_H
#define SIMDIALOG_H
//
#include "ui_simDlg.h"
#include <qwt_plot_curve.h>

//
class SimDialog : public QDialog, public Ui::Dialog
{
  Q_OBJECT
  public:
    SimDialog( QWidget * parent = 0, Qt::WFlags f = 0 );
    
  private slots:  
    void on_runButton_clicked();
    void on_saveButton_clicked();

  private:
    QwtPlotCurve *cReg;
    QwtPlotCurve *cOct;

};
#endif





