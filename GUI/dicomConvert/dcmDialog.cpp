#include "dcmDialog.h"
#include <QString>

#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"     
// #include "dcmDataSeries.h"
#include "saveVolumeDialog.h"

// analyze support ...
#include "dbh.h"

dcmDialog::dcmDialog(QWidget *parent)
: QDialog(parent)
{
  ui.setupUi(this);

  ui.m_pPreview->setColumnCount(6);
  // ui.m_pPreview->setIconSize(QSize(64, 64));
  // ui.m_pPreview->setIconSize(QSize(24, 24));

  QStringList headers;
  headers << tr("Subject") << tr("Protocol") << tr("Study ID") << tr("Series") << tr("Dimensions") << tr("kV");
  ui.m_pPreview->setHeaderLabels(headers);
  ui.m_pPreview->setSelectionMode(QAbstractItemView::ExtendedSelection);
	ui.m_pPreview->setSortingEnabled(true);
	ui.m_pPreview->sortByColumn(3);

  connect( ui.m_pPreview, SIGNAL(itemExpanded ( QTreeWidgetItem * )), this, SLOT(resizeCols()) );
  connect( ui.m_pPreview, SIGNAL(itemSelectionChanged ()), this, SLOT(updatePreview()) );

  connect( ui.m_pDirBrowseBut, SIGNAL( clicked() ), this, SLOT( changeDir()));
  // connect( ui.m_pScanBut, SIGNAL( clicked() ), this, SLOT( scanDir()));
  connect( ui.m_pScanBut, SIGNAL( clicked() ), this, SLOT( scanDirRecursive()));
  connect( ui.m_pSave3dBut, SIGNAL( clicked() ), this, SLOT( save3D()));
  connect( ui.m_pSave4dBut, SIGNAL( clicked() ), this, SLOT( save4D()));
}

void dcmDialog::changeDir() {
  // get new directory ... and set it in the list view ...
  QString s = QFileDialog::getExistingDirectory(
                    this,
                    "Choose a directory",
                    QDir::currentPath(),
                    QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  // now into the comboBox.
  ui.m_pDcmDir->addItem(s);
  ui.m_pDcmDir->setCurrentIndex(ui.m_pDcmDir->count() - 1);

}

bool dcmDialog::scanDirRecursive() {
  
  // re-implemented from dcmDialog.

  // need a list of stringlists ...
  QStringList patients;
  QStringList files;

  dseries.clear();  
  QString patient, protocol, ppos; 
  int study, series;
  int x,y;

  int numFiles = 0, i=0;

  QDir directory( ui.m_pDcmDir->currentText() );

  // need to scan the directory recursively and get the total number and list of files ...
  scanDir(directory.absolutePath(), files);
  numFiles = files.size();

  // qDebug("HS: Files in all subdirs: %d", numFiles);

  // std::ofstream out("files.txt");

  QProgressDialog progress("Scanning DICOM directory...", "Abort Scan", 0, numFiles, this);
  progress.setValue(i);
  foreach(QString fname, files) {
    progress.setValue(i++);
    // progress.setLabelText(filestring);
    qApp->processEvents();

    if (progress.wasCanceled())
      break;
        
    // filenames to avoid (mostly directories) ...
    if ( fname.contains(".calib", Qt::CaseInsensitive) || fname.contains(".avi", Qt::CaseInsensitive) || fname.contains("thumbnails", Qt::CaseInsensitive) ||  fname.contains("dti", Qt::CaseInsensitive) ||  fname.contains("dicomdir", Qt::CaseInsensitive) ) 
  
      continue;

    // out << qPrintable(fname) << std::endl;
    // std::cout << qPrintable(fname) << std::endl;
    // continue;

    DcmFileFormat fileformat;
    OFCondition status = fileformat.loadFile(fname.toLatin1());
    if (status.good())
    {
      OFString tmpStr;
      if (fileformat.getDataset()->findAndGetOFString(DCM_PatientsName, tmpStr).good()) {
        patient = tmpStr.c_str ();
        // clean up up the patient name ...
        patient = patient.toLower();
        patient = patient.split("^").join(", ");
        // patient.replace(QRegExp("^(.*), (.*)$"), "\\2 \\1");

        if ( ! patients.contains(patient, Qt::CaseInsensitive) )
          patients.append(patient);

      } else {
				patient = "MissingTag";
				if ( ! patients.contains(patient, Qt::CaseInsensitive) )
          patients.append(patient);
        qDebug("Error: cannot access Patient's Name!");
      }
      if (fileformat.getDataset()->findAndGetOFString(DCM_ProtocolName, tmpStr).good()) {
        protocol = tmpStr.c_str ();
      } else {
        protocol = "MissingTag";
				qDebug("Error: cannot access Protocol Name!");
      }
      if (fileformat.getDataset()->findAndGetOFString(DCM_StudyID, tmpStr).good()) {
        study = atoi(tmpStr.c_str ());
      } else {
				study = 1;
        qDebug("Error: cannot access Study Name!");
      } 
      if (fileformat.getDataset()->findAndGetOFString(DCM_SeriesNumber, tmpStr).good()) {
        series = atoi(tmpStr.c_str ());
      } else {
				series = 1;
        qDebug("Error: cannot access Series Name!");
      }

	  if (fileformat.getDataset()->findAndGetOFString(DCM_SliceLocation, tmpStr).good()) {
        ppos = QString(tmpStr.c_str ());
      } else {
				ppos = "0";
        qDebug("Error: cannot access Patient position!");
      }

		// kV
		unsigned int kvp;

		if (fileformat.getDataset()->findAndGetOFString(DCM_KVP, tmpStr).good()) {
			kvp = atoi(tmpStr.c_str ());
		} else {
			kvp = 0;
			qDebug("Error: cannot access XRay KVP !");
		}


      dcmDataSeries dsr(patient, protocol, study, series, ppos, kvp);
      // search to see if dsr is already in the list ...
      int idx = dseries.indexOf(dsr); 
      if ( idx == -1){
        dsr.files << fname;
        // DIMENSIONS
        if (fileformat.getDataset()->findAndGetOFString(DCM_Rows, tmpStr).good()) {
          y = atoi(tmpStr.c_str ());
        } else {
					y=0;
          qDebug("Error: cannot access Image rows !");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_Columns, tmpStr).good()) {
          x = atoi(tmpStr.c_str ());
        } else {
					x=0;
          qDebug("Error: cannot access Image Columns!");
        }

        dsr.x = x; dsr.y = y;

				// data type ...
				unsigned int ba;

				if (fileformat.getDataset()->findAndGetOFString(DCM_BitsAllocated, tmpStr).good()) {
					ba = atoi(tmpStr.c_str ());
				} else {
					ba=8;
					qDebug("Error: cannot access Image bits allocated !");
				}

        // preview ...
        unsigned char *prvu = new unsigned char[4*x*y];
        
				if (ba > 8) {
					const Uint16 *pixelData = NULL;
					unsigned long count = 0;
					if (fileformat.getDataset()->findAndGetUint16Array(DCM_PixelData, pixelData, &count).good())
					{ 
						int _max = 0;
						for (int i=0; i<x*y; i++)
							if (_max < pixelData[i]) 
								_max = pixelData[i];
						float fac = 255.0/_max;

						// std::cout << "Max is " << _max << " and fac is " << fac << std::endl;

						/* convert the pixel data */

						for (int i=0; i<x*y; i++) {
							prvu[4*i] = (unsigned char ) (fac*pixelData[i] ) ;
							prvu[4*i+1] = (unsigned char) (fac*pixelData[i]) ;
							prvu[4*i+2] = (unsigned char) (fac*pixelData[i]) ;
							prvu[4*i+3] = (unsigned char) (fac*pixelData[i]) ;
						}
					} else {
						qDebug("Error reading pixel data");
					}
				} else {
					const Uint8 *pixelData = NULL;
					unsigned long count = 0;
					if (fileformat.getDataset()->findAndGetUint8Array(DCM_PixelData, pixelData, &count).good())
					{ 
						for (int i=0; i<x*y; i++) {
							prvu[4*i] = pixelData[i]  ;
							prvu[4*i+1] = pixelData[i] ;
							prvu[4*i+2] = pixelData[i] ;
							prvu[4*i+3] = pixelData[i] ;
						}
					} else {
						qDebug("Error reading pixel data");
					}
				}

        // set this 
        dsr.preview = prvu;

        dseries.append(dsr);
      } else { 
        dseries[idx].files << fname;
      }
    } else
      qDebug("Error: cannot read DICOM file (%s)", status.text());
  }
  progress.setValue(numFiles);
  // out.close();

  // instead of std::out .. insert them into the list view ...

  ui.m_pPreview->clear();
 
  QList<QTreeWidgetItem*> _patientList;

  // insert all patients as top level widgetitems ...
  for (i=0; i<patients.size(); i++) {
	QTreeWidgetItem *pat = new QTreeWidgetItem( ui.m_pPreview );
    pat->setText(0, patients[i]);
    pat->setIcon(0, QIcon(QPixmap(":/res/images/patient.png")));
    // pat->setText(1, QString("dummy"));
    _patientList.append(pat);
  }

  // QTreeWidgetItem *lastItem;

  // now enter all the protocols ...
	char my_series[100];
  for (i=0; i<dseries.size(); i++) {
	  int _pid = patients.indexOf(dseries[i].patient);
	  QTreeWidgetItem *srs = new QTreeWidgetItem( _patientList[_pid] );
	  srs->setText(0, dseries[i].patient);
	  srs->setText(1, dseries[i].protocol);
	  srs->setText(2, QString::number(dseries[i].study));
	  sprintf(my_series, "%.2d", dseries[i].series);
		srs->setText(3, QString(my_series));
		// srs->setText(3, QString::number(dseries[i].series));
	  srs->setText(4, QString("%1x%2x%3").arg(dseries[i].x).arg(dseries[i].y).arg(dseries[i].files.size()));
		srs->setText(5, QString::number(dseries[i].kvp));
	  // srs->setText(5, QString::number(frameCnt));
	  // show a preview as the icon ...
	  QImage prv( dseries[i].preview, dseries[i].x, dseries[i].y, QImage::Format_RGB32 );
	  srs->setData(0, Qt::ToolTipRole, QPixmap::fromImage(prv) );
	  srs->setData(0, Qt::WhatsThisRole, dseries[i].ppos );

	  // free memory ...
	  // delete [] dseries[i].preview;

	  srs->setIcon(0, QIcon(QPixmap(":/res/images/data.png"))); 
	  // lastItem = srs;
	  // srs->setIcon(0, QIcon(QPixmap::fromImage(prv ))); 
	  // ui.m_lblPreview->setPixmap( QPixmap::fromImage( prv ) );
  }
  return true;
}

void dcmDialog::scanDir() {
  // get the current Index ... 
  dcmDir.setPath( ui.m_pDcmDir->currentText() );
  QStringList qsl;
  qsl += "*";
  // need a list of stringlists ...
  QStringList patients;

  dseries.clear();  

  QString fname;
  QString patient, protocol, ppos; 
  int study, series;
  int x,y;

  int numFiles = 0, i=0;

  foreach(QString filestring, dcmDir.entryList(qsl, QDir::Files))
    numFiles++;

  qDebug("Files in dir: %d", numFiles);

  QProgressDialog progress("Scanning DICOM directory...", "Abort Scan", 0, numFiles, this);
  progress.setValue(i);
  foreach(QString filestring, dcmDir.entryList(qsl, QDir::Files)) {
      progress.setValue(i++);
      // progress.setLabelText(filestring);
      fname = dcmDir.filePath(filestring);
      qApp->processEvents();

      if (progress.wasCanceled())
        break;

      DcmFileFormat fileformat;
      OFCondition status = fileformat.loadFile(fname.toLatin1());
      if (status.good())
      {
        OFString tmpStr;
        if (fileformat.getDataset()->findAndGetOFString(DCM_PatientsName, tmpStr).good()) {
          patient = tmpStr.c_str ();
          // clean up up the patient name ...
          patient = patient.toLower();
          patient = patient.split("^").join(", ");
          // patient.replace(QRegExp("^(.*), (.*)$"), "\\2 \\1");
          
          if ( ! patients.contains(patient, Qt::CaseInsensitive) )
            patients.append(patient);

        } else {
          qDebug("Error: cannot access Patient's Name!");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_ProtocolName, tmpStr).good()) {
          protocol = tmpStr.c_str ();
        } else {
          qDebug("Error: cannot access Protocol Name!");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_StudyID, tmpStr).good()) {
          study = atoi(tmpStr.c_str ());
        } else {
          qDebug("Error: cannot access Study Name!");
        } 
        if (fileformat.getDataset()->findAndGetOFString(DCM_SeriesNumber, tmpStr).good()) {
          series = atoi(tmpStr.c_str ());
        } else {
          qDebug("Error: cannot access Series Name!");
        }
		if (fileformat.getDataset()->findAndGetOFString(DCM_SliceLocation, tmpStr).good()) {
			ppos = QString(tmpStr.c_str ());
		} else {
			qDebug("Error: cannot access Patient position!");
		}
        
		// kV
		unsigned int kvp;

		if (fileformat.getDataset()->findAndGetOFString(DCM_KVP, tmpStr).good()) {
			kvp = atoi(tmpStr.c_str ());
		} else {
			kvp = 0;
			qDebug("Error: cannot access XRay KVP !");
		}

        dcmDataSeries dsr(patient, protocol, study, series, ppos, kvp);
        // search to see if dsr is already in the list ...
        int idx = dseries.indexOf(dsr); 
        if ( idx == -1){
          dsr.files << filestring;
          // DIMENSIONS
          if (fileformat.getDataset()->findAndGetOFString(DCM_Rows, tmpStr).good()) {
            y = atoi(tmpStr.c_str ());
          } else {
            qDebug("Error: cannot access Image rows !");
          }
          if (fileformat.getDataset()->findAndGetOFString(DCM_Columns, tmpStr).good()) {
            x = atoi(tmpStr.c_str ());
          } else {
            qDebug("Error: cannot access Image Columns!");
          }

          dsr.x = x; dsr.y = y;

					// data type ...
					unsigned int ba;

					if (fileformat.getDataset()->findAndGetOFString(DCM_BitsAllocated, tmpStr).good()) {
						ba = atoi(tmpStr.c_str ());
					} else {
						ba=8;
						qDebug("Error: cannot access Image bits allocated !");
					}

					// preview ...
					unsigned char *prvu = new unsigned char[4*x*y];

					if (ba > 8) {
						const Uint16 *pixelData = NULL;
						unsigned long count = 0;
						if (fileformat.getDataset()->findAndGetUint16Array(DCM_PixelData, pixelData, &count).good())
						{ 
							int _max = 0;
							for (int i=0; i<x*y; i++)
								if (_max < pixelData[i]) 
									_max = pixelData[i];
							float fac = 255.0/_max;

							// std::cout << "Max is " << _max << " and fac is " << fac << std::endl;

							/* convert the pixel data */

							for (int i=0; i<x*y; i++) {
								prvu[4*i] = (unsigned char ) (fac*pixelData[i] ) ;
								prvu[4*i+1] = (unsigned char) (fac*pixelData[i]) ;
								prvu[4*i+2] = (unsigned char) (fac*pixelData[i]) ;
								prvu[4*i+3] = (unsigned char) (fac*pixelData[i]) ;
							}
						} else {
							qDebug("Error reading pixel data");
						}
					} else {
						const Uint8 *pixelData = NULL;
						unsigned long count = 0;
						if (fileformat.getDataset()->findAndGetUint8Array(DCM_PixelData, pixelData, &count).good())
						{ 
							for (int i=0; i<x*y; i++) {
								prvu[4*i] = pixelData[i]  ;
								prvu[4*i+1] = pixelData[i] ;
								prvu[4*i+2] = pixelData[i] ;
								prvu[4*i+3] = pixelData[i] ;
							}
						} else {
							qDebug("Error reading pixel data");
						}
					}

          // set this 
          dsr.preview = prvu;

          dseries.append(dsr);
        } else { 
          dseries[idx].files << filestring;
        }
      } else
        qDebug("Error: cannot read DICOM file (%s)", status.text());
  }
  progress.setValue(numFiles);

  // instead of std::out .. insert them into the list view ...
  ui.m_pPreview->clear();
  QList<QTreeWidgetItem*> _patientList;
    
  // insert all patients as top level widgetitems ...
  for (i=0; i<patients.size(); i++) {
    QTreeWidgetItem *pat = new QTreeWidgetItem(ui.m_pPreview);
    pat->setText(0, patients[i]);
    pat->setIcon(0, QIcon(QPixmap(":/res/images/patient.png")));
    // pat->setText(1, QString("dummy"));
    _patientList.append(pat);
  }

  // now enter all the protocols ...
	char my_series[100];
  for (i=0; i<dseries.size(); i++) {
    int _pid = patients.indexOf(dseries[i].patient);
    QTreeWidgetItem *srs = new QTreeWidgetItem( _patientList[_pid] );
    srs->setText(0, dseries[i].patient);
    srs->setText(1, dseries[i].protocol);
    srs->setText(2, QString::number(dseries[i].study));
    sprintf(my_series, "%.2d", dseries[i].series);
		srs->setText(3, QString(my_series));
		// srs->setText(3, QString::number(dseries[i].series));
    srs->setText(4, QString::number(dseries[i].files.size()));
		srs->setText(5, QString::number(dseries[i].kvp));
    // show a preview as the icon ...
    QImage prv( dseries[i].preview, dseries[i].x, dseries[i].y, QImage::Format_RGB32 );
    srs->setData(0, Qt::ToolTipRole, QPixmap::fromImage(prv) );
    
    srs->setIcon(0, QIcon(QPixmap(":/res/images/data.png"))); 
    // srs->setIcon(0, QIcon(QPixmap::fromImage(prv ))); 
    // ui.m_lblPreview->setPixmap( QPixmap::fromImage( prv ) );
  }

  // resize columns ...
  ui.m_pPreview->resizeColumnToContents(2);
  ui.m_pPreview->resizeColumnToContents(0);
  ui.m_pPreview->resizeColumnToContents(1);
}

void dcmDialog::resizeCols() {
  for (int i=0; i<5; i++)
    ui.m_pPreview->resizeColumnToContents(i);
}

void dcmDialog::save3D() {
  // Lets save 3D ... :) 

  QList<dcmDataSeries> _saveseries;
  // QList<dcmDataSeries> _series;
  dcmDataSeries dsr;
  bool ok;
  FILE *fp1;

  // first get the selection from m_pPreview ...
  QList<QTreeWidgetItem *> items = ui.m_pPreview->selectedItems();
  if ( ! items.size() ) {
    qDebug("No Items selected");
    return;
  }

  // this can be changed into a list of dcmDataSeries ...
  for (int i=0; i<items.size(); i++) {
    dsr.patient = items[i]->text(0);  
    dsr.protocol = items[i]->text(1);
    dsr.study = items[i]->text(2).toInt(&ok);  
    dsr.series = items[i]->text(3).toInt(&ok);
	QVariant ppos = items[i]->data(0, Qt::WhatsThisRole);
	dsr.ppos = ppos.value<QString>();
    int idx = dseries.indexOf(dsr);
    _saveseries.append(dseries[idx]);
  }

  // Get path to directory ... + prefix ... makes sense ?

  saveVolumeDialog dlg;
  int stat = dlg.exec();
  if (stat) {
    // now use dcmtk to stack the correct volumes together ...
    // std::cout << "Yes" << std::endl;
    // collect all the info required for the entire set first ...

    int x, y,z; // x/y/z dims ;
    double sx=1.0, sy=1.0, sz=1.0; // pixel spacing in the x,y and x dimensions ...
    int bitsAlloc, bitsStored;

    // @todo need to get offset and orientation too ...
   
    int cnt = _saveseries.size();
    QProgressDialog progress("Saving 3D volume(S)...", "Abort Scan", 0, cnt, this);
    progress.setValue(0);

    // std::cout << "Starting Loop" << std::endl;
    // outer loop over all subjects ...
    for (int t=0; t<cnt; t++) {
      progress.setLabelText("Saving Volume Frame "+ QString::number(t));

      z = _saveseries[t].files.size();
      // use the first slice of the frame to get most of the info ...
      QString fname = dcmDir.filePath(_saveseries[t].files[0]);
      DcmFileFormat fileformat;
      OFCondition status = fileformat.loadFile(fname.toLatin1());
      if (status.good())
      {
        OFString tmpStr;

        // First check if Transfer Syntax is all right ... 
        /*
           if (fileformat.getDataset()->findAndGetOFString(DCM_TransferSyntaxUID, tmpStr).good()) {
           std::cout << "Transfer Syntax: " << tmpStr;
           } else {
           qDebug("Error: cannot access Transfer Syntax !");
           } */

        // DIMENSIONS
        if (fileformat.getDataset()->findAndGetOFString(DCM_Rows, tmpStr).good()) {
          y = atoi(tmpStr.c_str ());
          std::cout << "Y " << tmpStr;
        } else {
          qDebug("Error: cannot access Image rows !");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_Columns, tmpStr).good()) {
          x = atoi(tmpStr.c_str ());
          std::cout << " X " << tmpStr;
        } else {
          qDebug("Error: cannot access Image Columns!");
        }

        // SPACING
        if (fileformat.getDataset()->findAndGetOFString(DCM_PixelSpacing, tmpStr).good()) {
          std::cout << " Spacing is " << tmpStr << std::endl;
          sx = sy = atof(tmpStr.c_str() );
        } else {
          qDebug("Error: cannot access Pixel Spacing!");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_SliceThickness, tmpStr).good()) {
          sz = atof(tmpStr.c_str() );
        } else {
          qDebug("Error: cannot access Pixel Spacing - (slice thickness)!");
        }

        // pixel info ...
        if (fileformat.getDataset()->findAndGetOFString(DCM_BitsAllocated, tmpStr).good()) {
          bitsAlloc = atoi(tmpStr.c_str() );
        } else {
          qDebug("Error: cannot access Pixel BitsAllocated!");
        }
        if (fileformat.getDataset()->findAndGetOFString(DCM_BitsStored, tmpStr).good()) {
          bitsStored = atoi(tmpStr.c_str() );
        } else {
          qDebug("Error: cannot access Pixel BitsStored!");
        }
      }

      // std::cout << "Bits: " << bitsAlloc << " " << bitsStored << std::endl;

      std::cout << "Got INFO" << std::endl;
      // Prepare the header based on the format ...
      if (dlg.getSaveFormat() == QString("MetaIO"))
        std::cout << "Write MHD" << std::endl;
      else 
        std::cout << "Write Analyze" << std::endl;

      // create the volInfo ...
      volInfo info;
      info.x = x; info.y = y; info.z = z;
      info.sx = sx; info.sy = sy; info.sz = sz;
      info.bs = 8; info.ba = 8;

      // @bug @todo for now I shall convert all unsigned short images to
      // normalized unsigned char images. I think works better this way.

      // allocate memory for a single volume ...
      unsigned char *vol = new unsigned char[x*y*z];
      unsigned short *vol2 = new unsigned short[x*y*z];

      saveHeader(VOL_FMT_ANALYZE, dlg.getSavePath() +"/" + dlg.getFilePrefix() + QString().sprintf(".%02d", t) + ".hdr", info);

      for (int k=0; k<z; k++) {
        // load the correct dataset ...
        fname = dcmDir.filePath(_saveseries[t].files[k]);
        // fname = dcmDir.filePath(_series[t].files[k]);
        status = fileformat.loadFile(fname.toLatin1());
        progress.setValue(t);
        qApp->processEvents();
        // std::cout << qPrintable(fname) << std::endl;

        // read in the slices ... 
        const Uint16 *pixelData = NULL;
        unsigned long count = 0;
        if (fileformat.getDataset()->findAndGetUint16Array(DCM_PixelData, pixelData, &count).good())
        {
          /* convert the pixel data */
          // std::cout << "Converting pizedl data" << std::endl;
          for (int i=0; i<x*y; i++)
            vol2[k*x*y + i] = pixelData[i];
        } else {
          qDebug("Error reading pixel data");
        }

      }
      // convert image to unsigned char ...  
      unsigned int _max =0;  
      for (int i=0; i<x*y*z; i++) {
        if ( _max < vol2[i])
           _max = vol2[i];
      }
      for (int i=0; i<x*y*z; i++) 
        vol[i] = static_cast<unsigned char>((250.0/_max)* vol2[i]);
      // write out image ...
      QString volname = dlg.getSavePath() +"/" + dlg.getFilePrefix() + QString().sprintf(".%02d", t) + ".img";
      fp1 = fopen( volname.toLatin1(), "wb");
      fwrite(vol, 1, x*y*z, fp1);
      fclose(fp1);
      delete [] vol;
      delete [] vol2;
    }
    progress.setValue(cnt);

  } // if stat ... 

}

void dcmDialog::save4D() {
  // Lets save 4D ... :) 

  QList<dcmDataSeries> _saveseries;
  // QList<dcmDataSeries> _series;
  dcmDataSeries dsr;
  bool ok;
  FILE *fp1;

  // first get the selection from m_pPreview ...
  QList<QTreeWidgetItem *> items = ui.m_pPreview->selectedItems();
  if ( ! items.size() ) {
    qDebug("No Items selected");
    return;
  }

  // this can be changed into a list of dcmDataSeries ...
  for (int i=0; i<items.size(); i++) {
    dsr.patient = items[i]->text(0);  
    dsr.protocol = items[i]->text(1);
    dsr.study = items[i]->text(2).toInt(&ok);  
    dsr.series = items[i]->text(3).toInt(&ok);
	QVariant ppos = items[i]->data(0, Qt::WhatsThisRole);
	dsr.ppos = ppos.value<QString>();
    int idx = dseries.indexOf(dsr);
    _saveseries.append(dseries[idx]);
  }

  // check to see if all of them have the same length ...
  int len = _saveseries[0].files.size();
  int z = _saveseries.size();

  for (int i=1; i<_saveseries.size(); i++) {
    if ( len != _saveseries[i].files.size() ) {
      // error message ...
      QMessageBox::critical(0, "Dicom Import",
          QString("The selected sequences do not have the same number of frames.\n\n") +
          QString("Please retry\n"));
      return;
    }
  }

  // Get path to directory ... + prefix ... makes sense ?

  saveVolumeDialog dlg;
  int stat = dlg.exec();
  if (stat) {
    // now use dcmtk to stack the correct volumes together ...
    std::cout << "Yes" << std::endl;
    // collect all the info required for the entire set first ...

    int x, y; // x/y dims ... we already know z and t;
    double sx=1.0, sy=1.0, sz=1.0; // pixel spacing in the x,y and x dimensions ...
    int bitsAlloc, bitsStored;

    // @todo need to get offset and orientation too ...

    // use the first slice of the first frame to get most of the info ...
    QString fname = dcmDir.filePath(_saveseries[0].files[0]);
    DcmFileFormat fileformat;
    OFCondition status = fileformat.loadFile(fname.toLatin1());
    if (status.good())
    {
      OFString tmpStr;

      // First check if Transfer Syntax is all right ... 
      /*
         if (fileformat.getDataset()->findAndGetOFString(DCM_TransferSyntaxUID, tmpStr).good()) {
         std::cout << "Transfer Syntax: " << tmpStr;
         } else {
         qDebug("Error: cannot access Transfer Syntax !");
         } */

      // DIMENSIONS
      if (fileformat.getDataset()->findAndGetOFString(DCM_Rows, tmpStr).good()) {
        y = atoi(tmpStr.c_str ());
        std::cout << "Y " << tmpStr;
      } else {
        qDebug("Error: cannot access Image rows !");
      }
      if (fileformat.getDataset()->findAndGetOFString(DCM_Columns, tmpStr).good()) {
        x = atoi(tmpStr.c_str ());
        std::cout << " X " << tmpStr;
      } else {
        qDebug("Error: cannot access Image Columns!");
      }

      // SPACING
      if (fileformat.getDataset()->findAndGetOFString(DCM_PixelSpacing, tmpStr).good()) {
        std::cout << " Spacing is " << tmpStr << std::endl;
        sx = sy = atof(tmpStr.c_str() );
      } else {
        qDebug("Error: cannot access Pixel Spacing!");
      }
      if (fileformat.getDataset()->findAndGetOFString(DCM_SliceThickness, tmpStr).good()) {
        sz = atof(tmpStr.c_str() );
      } else {
        qDebug("Error: cannot access Pixel Spacing - (slice thickness)!");
      }

      // pixel info ...
      if (fileformat.getDataset()->findAndGetOFString(DCM_BitsAllocated, tmpStr).good()) {
        bitsAlloc = atoi(tmpStr.c_str() );
      } else {
        qDebug("Error: cannot access Pixel BitsAllocated!");
      }
      if (fileformat.getDataset()->findAndGetOFString(DCM_BitsStored, tmpStr).good()) {
        bitsStored = atoi(tmpStr.c_str() );
      } else {
        qDebug("Error: cannot access Pixel BitsStored!");
      }
    }

    std::cout << "Bits: " << bitsAlloc << " " << bitsStored << std::endl;

    // Prepare the header based on the format ...
    if (dlg.getSaveFormat() == QString("MetaIO"))
      std::cout << "Write MHD" << std::endl;
    else 
      std::cout << "Write Analyze" << std::endl;

    // create the volInfo ...
    volInfo info;
    info.x = x; info.y = y; info.z = z;
    info.sx = sx; info.sy = sy; info.sz = sz;
    info.bs = 8; info.ba = 8;

    // @bug @todo for now I shall convert all unsigned short images to
    // normalized unsigned char images. I think works better this way.

    // allocate memory for a single volume ...
    unsigned char *vol = new unsigned char[x*y*z];

    /// compute the total number of images to be read so that the progress
    //counter can be set ... 
    int cnt = z*len, ii=0;

    QProgressDialog progress("Saving 4D volume...", "Abort Scan", 0, cnt, this);
    progress.setValue(ii);

    // outer loop over len ...
    for (int t=0; t<len; t++) {
      progress.setLabelText("Saving Volume Frame "+ QString::number(t));
      saveHeader(VOL_FMT_ANALYZE, dlg.getSavePath() +"/" + dlg.getFilePrefix() + QString().sprintf(".%02d", t) + ".hdr", info);
      for (int k=0; k<z; k++) {
        // load the correct dataset ...
        fname = dcmDir.filePath(_saveseries[k].files[t]);
        // fname = dcmDir.filePath(_series[t].files[k]);
        status = fileformat.loadFile(fname.toLatin1());
        progress.setValue(ii++);
        qApp->processEvents();
        std::cout << qPrintable(fname) << std::endl;
        
        // read in the slices ... 
        const Uint16 *pixelData = NULL;
        unsigned long count = 0;
        if (fileformat.getDataset()->findAndGetUint16Array(DCM_PixelData, pixelData, &count).good())
        {
          /* convert the pixel data */
          // std::cout << "Converting pizedl data" << std::endl;
          for (int i=0; i<x*y; i++)
            vol[k*x*y + i] = (unsigned char )( pixelData[i] >> 1);
        } else {
          qDebug("Error reading pixel data");
        }

      }
      // write out image ...
      QString volname = dlg.getSavePath() +"/" + dlg.getFilePrefix() + QString().sprintf(".%02d", t) + ".img";
      fp1 = fopen( volname.toLatin1(), "wb");
      fwrite(vol, 1, x*y*z, fp1);
      fclose(fp1);
    }
    progress.setValue(cnt);

    delete [] vol;
  } // if stat ... 
}

void dcmDialog::saveHeader(int format, QString fname, volInfo _info) {

  if (format == VOL_FMT_MHD ) {
	// Write out MHD header.	

  } else if (format == VOL_FMT_ANALYZE) {
    struct dsr hdr;
    FILE *fp;

    static char DataTypes[9][12] = {"UNKNOWN", "BINARY", "CHAR", "SHORT", "INT","FLOAT", "COMPLEX", "DOUBLE", "RGB"};
    static int DataTypeSizes[9] = {0,1,8,16,32,32,64,64,24};

    // std::cout << "Saving to file : " << qPrintable(fname) << std::endl;

    memset(&hdr,0, sizeof(struct dsr));
    for(int i=0;i<8;i++)
      hdr.dime.pixdim[i] = 0.0;

    hdr.dime.vox_offset = 0.0;
    hdr.dime.funused1 = 0.0;
    hdr.dime.funused2 = 0.0;
    hdr.dime.funused3 = 0.0;
    hdr.dime.cal_max = 0.0;
    hdr.dime.cal_min = 0.0;
    hdr.dime.datatype = -1;

    for(int i=1;i<=8;i++)
      if(!strcmp("CHAR", DataTypes[i]))
      {
        hdr.dime.datatype = (1<<(i-1));
        hdr.dime.bitpix = DataTypeSizes[i];
        break;
      }
    if((fp=fopen(fname.toLatin1(),"w"))==0)
    {
      qDebug("unable to create header file");
      return;
    }
    hdr.dime.dim[0] = 4; /* all Analyze images are taken as 4 dimensional */
    hdr.hk.regular = 'r';
    hdr.hk.sizeof_hdr = sizeof(struct dsr);
    hdr.dime.dim[1] = _info.x; /* slice width in pixels */
    hdr.dime.dim[2] = _info.y; /* slice height in pixels */
    hdr.dime.dim[3] = _info.z; /* volume depth in slices */
    hdr.dime.dim[4] = 1; /* number of volumes per file */
    hdr.dime.glmax = 0; /* maximum voxel value */
    hdr.dime.glmin = 255; /* minimum voxel value */

    /* Set the voxel dimension fields:
       A value of 0.0 for these fields implies that the value is unknown.
       Change these values to what is appropriate for your data
       or pass additional command line arguments */
    hdr.dime.pixdim[1] = _info.sx; /* voxel x dimension */
    hdr.dime.pixdim[2] = _info.sy; /* voxel y dimension */
    hdr.dime.pixdim[3] = _info.sz; /* pixel z dimension, slice thickness */

    /* Assume zero offset in .img file, byte at which pixel
       data starts in the image file */
    hdr.dime.vox_offset = 0.0;
    /* Planar Orientation; */
    /* Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
       Movie flag ON: 3 = transverse, 4 = coronal, 5 = sagittal */
    hdr.hist.orient = 0;
    /* up to 3 characters for the voxels units label; i.e. mm., um., cm. */ 
    strcpy(hdr.dime.vox_units," ");
    /* up to 7 characters for the calibration units label; i.e. HU */
    strcpy(hdr.dime.cal_units," ");
    /* Calibration maximum and minimum values;
       values of 0.0 for both fields imply that no
       calibration max and min values are used */
    hdr.dime.cal_max = 0.0;
    hdr.dime.cal_min = 0.0;
    fwrite(&hdr,sizeof(struct dsr),1,fp);
    fclose(fp);
  } else {
    qDebug("Unknown output format");
  }
}

void dcmDialog::updatePreview() {
  QList<QTreeWidgetItem *> items = ui.m_pPreview->selectedItems();
  // only update on single item ...
  if (items.size() == 1) {
    QVariant pix = items[0]->data(0, Qt::ToolTipRole);
    // update the preview ..
    ui.m_lblPreview->setPixmap(pix.value<QPixmap>());
  }
}

// Private utility functions ...
void dcmDialog::addFilesFromDir( const QString& dcmDir, QStringList& fileList) {
  QDir theDir( dcmDir );
  QStringList localList;

  theDir.setFilter(QDir::Files);
  foreach(QString filestring, theDir.entryList()) {
    // skip DICOMDIR for now ....
    if (filestring.toLower() == QString("dicomdir"))
      continue;	
    fileList << theDir.filePath(filestring);
  }
}

void dcmDialog::scanDir(const QString& dcmDir, QStringList& fileList) {
  QString currentPath;
  QStringList dirList;
  QDir dir(dcmDir);

  // Stop recursion if the directory doesn't exist.
  if (!dir.exists()) {
    return;
  }

  dir.setSorting(QDir::Name);
  dir.setFilter(QDir::Dirs);
  dirList = dir.entryList();

  // Remove '.' and  '..'
  dirList.erase(dirList.begin());
  dirList.erase(dirList.begin());

  // Recurse through all directories
  QStringList::Iterator it;
  for( it = dirList.begin(); it != dirList.end(); ++it)
  {
    // Recursive call to fetch subdirectories
    currentPath = dcmDir + "//" + (*it);
    scanDir( currentPath, fileList );    
  } 

  // Add files in THIS directory as well.
  addFilesFromDir( dcmDir, fileList );
}