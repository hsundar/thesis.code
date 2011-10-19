#include <stdlib.h>
#include <qwt_painter.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_scale_widget.h>
#include <qwt_legend.h>
#include <qwt_scale_draw.h>
#include <qwt_math.h>
#include "data_plot.h"

#include <iostream>
//
//  Initialize main window
//
DataPlot::DataPlot(QWidget *parent):
    QwtPlot(parent)
{
    // Disable polygon clipping
    QwtPainter::setDeviceClipping(false);

    // We don't need the cache here
    canvas()->setPaintAttribute(QwtPlotCanvas::PaintCached, false);
    canvas()->setPaintAttribute(QwtPlotCanvas::PaintPacked, false);

#if QT_VERSION >= 0x040000
#ifdef Q_WS_X11
    /*
       Qt::WA_PaintOnScreen is only supported for X11, but leads
       to substantial bugs with Qt 4.2.x/Windows
     */
    canvas()->setAttribute(Qt::WA_PaintOnScreen, true);
#endif
#endif

    alignScales();
    
    /*
    m_size = 41;
    d_x = new double[m_size];
    d_y = new double[m_size];
    d_z = new double[m_size];
    //  Initialize data
    for (int i = -20; i< 21; i++)
    {
        d_x[i+20] = 0.5 * i;     // time axis
        d_y[i+20] = sin(i);
        d_z[i+20] = cos(i);
    }
    */
    // Assign a title
    setTitle("Comparing Similarity Measures");
    insertLegend(new QwtLegend(), QwtPlot::BottomLegend);

    // Insert new curves
    cReg = new QwtPlotCurve("Regular Similarity Measure");
    cReg->attach(this);

    cOct = new QwtPlotCurve("Octree based Similarity Measure");
    cOct->attach(this);

    // Set curve styles
    cReg->setPen(QPen(Qt::red, 2));
    cOct->setPen(QPen(Qt::blue, 2));

    d_x = NULL;
    d_y = NULL;
    d_z = NULL;
    /*
    // Attach (don't copy) data. Both curves use the same x array.
    cRight->setRawData(d_x, d_y, m_size);
    cLeft->setRawData(d_x, d_z, m_size);

    */
#if 1
    //  Insert zero line at y = 0
    QwtPlotMarker *mY = new QwtPlotMarker();
    mY->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
    mY->setLineStyle(QwtPlotMarker::HLine);
    mY->setYValue(0.0);
    mY->attach(this);
    
    //  Insert zero line at x = 0
    QwtPlotMarker *mX = new QwtPlotMarker();
    mX->setLabelAlignment(Qt::AlignRight|Qt::AlignTop);
    mX->setLineStyle(QwtPlotMarker::VLine);
    mX->setXValue(0.0);
    mX->attach(this);
#endif

    // Axis 
    setAxisTitle(QwtPlot::xBottom, "Variations");
    setAxisScale(QwtPlot::xBottom, -10, 10);

    setAxisTitle(QwtPlot::yLeft, "Similarity");
    setAxisScale(QwtPlot::yLeft, -1.5, 1.5);
}

//
//  Set a plain canvas frame and align the scales to it
//
void DataPlot::alignScales()
{
    // The code below shows how to align the scales to
    // the canvas frame, but is also a good example demonstrating
    // why the spreaded API needs polishing.

    canvas()->setFrameStyle(QFrame::Box | QFrame::Plain );
    canvas()->setLineWidth(1);

    for ( int i = 0; i < QwtPlot::axisCnt; i++ )
    {
        QwtScaleWidget *scaleWidget = (QwtScaleWidget *)axisWidget(i);
        if ( scaleWidget )
            scaleWidget->setMargin(0);

        QwtScaleDraw *scaleDraw = (QwtScaleDraw *)axisScaleDraw(i);
        if ( scaleDraw )
            scaleDraw->enableComponent(QwtAbstractScaleDraw::Backbone, false);
    }
}

void DataPlot::setPlotRange(int range, double step, QString leg) {
  
  m_range = range;
  m_step  = step;
    
  int _min = -range*step;
  int _max = range*step;
  std::cout << "Allocating memory" << std::endl;

  unsigned int _size = 2*range+1;
  if ( d_x != NULL) {
    delete [] d_x;
    d_x = NULL;
  }
  if ( d_y != NULL) {
    delete [] d_y;
    d_y = NULL;
  }
  if ( d_z != NULL) {
    delete [] d_z;
    d_z = NULL;
  }
  d_x = new double [_size];
  d_y = new double [_size];
  d_z = new double [_size];

  for (unsigned int i=0; i < _size; i++) {
    d_x[i] = m_step*i + _min;
    d_z[i] = 0.;
    d_y[i] = 0.;

  }

  std::cout << "setting raw data ..." << std::endl;
  cReg->setRawData(d_x, d_y, _size);
  cOct->setRawData(d_x, d_z, _size);
  std::cout << "finished setting raw data ..." << std::endl;
  
  setAxisTitle(QwtPlot::xBottom, leg);
  setAxisScale(QwtPlot::xBottom, _min, _max);
  std::cout << "Replotting ..." << std::endl;
  replot();
}

void DataPlot::setValue(unsigned int i, double val, bool useOct) {
  if (useOct)
    d_z[i] = val;
  else
    d_y[i] = val;
  replot();
}

void DataPlot::setRange(int r) {
  m_range = r;
}

void DataPlot::setStep(double s) {
  m_step = s;

}
