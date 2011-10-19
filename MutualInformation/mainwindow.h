#ifndef __MAINWINDOW_H_
#define __MAINWINDOW_H_

#include <QtGui>
#include "data_plot.h"

#include "Volume.h"
#include "oct_sim.h"

class MainWindow: public QMainWindow
{
  Q_OBJECT

  public:
    MainWindow();

  public slots:
    void onBrowseTarget();
    void onBrowseSource();
    
    void onRun();
    void onSave();

  protected:
    // Gui Objects
    DataPlot *m_plot;

    QCheckBox *useOctCheck;
    QLineEdit *targetFileName;
    QLineEdit *sourceFileName;
    
    QSpinBox *rangeSpin;
    QDoubleSpinBox *stepSpin;
    
    QComboBox *metricCombo;
    QComboBox *modesCombo;

    QPushButton *runButton;
    QPushButton *saveButton;

    // Processing Objects
    Volume *sourceVol;
    Volume *targetVol;

    oct_sim *sim; 
};



#endif

