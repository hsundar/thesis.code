#ifndef _DATA_PLOT_H
#define _DATA_PLOT_H 1

#include <qwt_plot.h>
#include <qwt_plot_curve.h>

class DataPlot : public QwtPlot
{
  Q_OBJECT

  public:
    DataPlot(QWidget* = NULL);

    void setPlotRange(int range, double step, QString leg); 
    void setValue(unsigned int i, double val, bool useOct);

  public slots:
    void setRange(int r);
    void setStep(double s);

  private:
      void alignScales();
      
      QwtPlotCurve *cReg;
      QwtPlotCurve *cOct;
      
      double *d_x; 

      double *d_y; 
      double *d_z;

      int 		m_range;
      double       	m_step;
      QString		m_xLabel;

};

#endif
