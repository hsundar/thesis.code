#ifndef __DCM_DATA_SERIES_H_
#define __DCM_DATA_SERIES_H_

#include <QtCore>

class dcmDataSeries {
public:
  dcmDataSeries() {};
  dcmDataSeries(QString p, QString q, int st, int sr):
    patient(p),
    protocol(q),
    study(st),
    series(sr)
  {};

  ~dcmDataSeries() {};

  inline bool operator != (const dcmDataSeries &other) { return ( ( patient != other.patient ) 
                                                               || ( series  != other.series  ) 
                                                               || ( study   != other.study   ) ); };
  inline bool operator == (const dcmDataSeries &other) { return ( ( patient == other.patient ) 
                                                               && ( series  == other.series  ) 
                                                               && ( study   == other.study   ) ); };
  QStringList files;

// protected:
  QString patient;
  QString protocol;
  int study;
  int series;
  int x,y;
  unsigned char *preview;
};


#endif
