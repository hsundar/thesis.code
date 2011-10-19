#include "simDialog.h"

#include <QtCore>
#include <QtGui>

#include <qwt_legend.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_marker.h>

#include "mainwindow.h"
#include "viewer3d.h"
#include "Volume.h"

//
SimDialog::SimDialog( QWidget * parent, Qt::WFlags f) 
	: QDialog(parent, f)
{
  setupUi(this);

  // setup the plot ...
  m_plot->setTitle("Comparing Similarity Measures");
  m_plot->insertLegend(new QwtLegend(), QwtPlot::BottomLegend);

  // set up the axes ...
  m_plot->setAxisTitle(QwtPlot::xBottom, "Variations");
  m_plot->setAxisScale(QwtPlot::xBottom, -10, 10);

  m_plot->setAxisTitle(QwtPlot::yLeft, "Similarity");
  m_plot->setAxisScale(QwtPlot::yLeft, -1.5, 1.5);

  // enable grid ...
  QwtPlotGrid *grid = new QwtPlotGrid;
  grid->enableXMin(true);
  grid->enableYMin(true);
  grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
  grid->setMinPen(QPen(Qt::gray, 0 , Qt::DotLine));
  grid->attach(m_plot);

  //  Insert zero line at y = 0
  QwtPlotMarker *mY = new QwtPlotMarker();
  mY->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
  mY->setLineStyle(QwtPlotMarker::HLine);
  mY->setYValue(0.0);
  mY->attach(m_plot);

  //  Insert zero line at x = 0
  QwtPlotMarker *mX = new QwtPlotMarker();
  mX->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
  mX->setLineStyle(QwtPlotMarker::VLine);
  mX->setXValue(0.0);
  mX->attach(m_plot);

  // Insert new curves
  cReg = new QwtPlotCurve("Regular Similarity Measure");
  cReg->attach(m_plot);

  cOct = new QwtPlotCurve("Octree based Similarity Measure");
  cOct->attach(m_plot);

  // Set curve styles
  cReg->setPen(QPen(Qt::red, 2));
  cOct->setPen(QPen(Qt::blue, 2));

  // Populate the volume combos ...
  MainWindow* _win  = (MainWindow *)this->parent();
  viewer3d *  _view = (viewer3d *) _win->getViewer();
  QList<Volume*> _vols = _view->getVolumes();

  // append Volume names to combo boxes ...
  // Need to modify Volume class so that it stored patient name ...
  // More importantly modify PatientBrowser so that the data can be directly
  // loaded.
}
//

void SimDialog::on_runButton_clicked()
{
  // TODO
}

void SimDialog::on_saveButton_clicked()
{
#if QT_VERSION >= 0x040100

#ifndef QT_NO_FILEDIALOG
  const QString fileName = QFileDialog::getSaveFileName(
      this, "Save File Name", QString(),
      "PDF Documents (*.pdf)");
#else
  const QString fileName = "similarity.pdf";
#endif

  if ( !fileName.isEmpty() )
  {
    QPrinter printer;
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOrientation(QPrinter::Landscape);
    printer.setOutputFileName(fileName);

    printer.setCreator("Octree Similarity");
    m_plot->print(printer);
  }
#endif
}

