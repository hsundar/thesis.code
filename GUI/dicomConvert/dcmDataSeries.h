#ifndef __DCM_DATA_SERIES_H_
#define __DCM_DATA_SERIES_H_

#include <QtCore>

class dcmDataSeries {
public:
  dcmDataSeries() {};
  dcmDataSeries(QString p, QString q, int st, int sr, QString pp, int kv):
    patient(p),
    protocol(q),
    study(st),
    series(sr),
		ppos(pp),
		kvp(kv)
  {};

  ~dcmDataSeries() {};

	/* 
	 Special modification made for tagged acquisitions from Siemens which keeps different slices in the same series.
	 */


  inline bool operator != (const dcmDataSeries &other) { 
	  if ( protocol.contains( "grid", Qt::CaseInsensitive) || protocol.contains( "Missing", Qt::CaseInsensitive) ) {
		  return ( ( patient != other.patient ) 
			  || ( series  != other.series  ) 
			  || ( study   != other.study   ) 
			  || ( ppos   != other.ppos   ) );
	  } else {
		  return ( ( patient != other.patient ) 
			  || ( series  != other.series  ) 
			  || ( study   != other.study   ) );
	  }
  };

  inline bool operator == (const dcmDataSeries &other) { 
	  if ( protocol.contains( "grid", Qt::CaseInsensitive) || protocol.contains( "Missing", Qt::CaseInsensitive) ) {
		  return ( ( patient == other.patient ) 
			  && ( series  == other.series  ) 
			  && ( study   == other.study   ) 
			  && ( ppos   == other.ppos   ) );
	  } else {
		  return ( ( patient == other.patient ) 
			  && ( series  == other.series  ) 
			  && ( study   == other.study   ) );
	  }
  };

  QStringList files;

// protected:
  QString patient;
  QString protocol;
	int kvp;
  int study;
  int series;
  QString ppos;
  int x,y;
  unsigned char *preview;
};


#endif
