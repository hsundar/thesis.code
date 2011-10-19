#include "patientBrowser.h"
#include <QFile>
#include <QMessageBox>
#include <QProgressDialog>
#include <QStringList>

// DCM includes ...
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dctk.h"  

// #include "dcmDataSeries.h"
//
patientBrowser::patientBrowser( QWidget * parent, Qt::WFlags f) : QDialog(parent, f)
{
  setupUi(this);
  directory = QDir();

  connect( patientTreeWidget, SIGNAL(itemExpanded ( QTreeWidgetItem * )), this, SLOT(resizeCols()) );
  connect( patientTreeWidget, SIGNAL(itemSelectionChanged ()), this, SLOT(updatePreview()) );
}

bool patientBrowser::initDataDirectory(QDir dir) {
  // First check if directory exists ...
  if (! dir.exists() )
    return false;

  directory = dir;
  QDir tn(dir);
  tn.cd("thumbnails");
  patientTreeWidget->setThumbnailDirectory(tn);

  // see if database file can be found here ...
  if ( dir.exists("cardiac.xml") ) {
    QFile file(directory.absoluteFilePath("cardiac.xml"));
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::warning(this, tr("Patient Browser"),
          tr("Cannot read file cardiac.xml:\n%a.")
          .arg(file.errorString()));
      return false;
    }
    return patientTreeWidget->read(&file);
  } else {	
    // if not .. force a scan of the folders to populate the folder.
    rescanDataDirectory();
    // now write to a file ...
    QFile file(directory.absoluteFilePath("cardiac.xml"));
    file.open(QFile::WriteOnly | QFile::Text);
    return patientTreeWidget->write(&file);

  }

  patientTreeWidget->resizeColumnToContents(0);
  return true;
}
//
bool patientBrowser::rescanDataDirectory() {
  
  // re-implemented from dcmDialog.

  // need a list of stringlists ...
  QStringList patients;
  QStringList files;

  dseries.clear();  
  QString patient, protocol; 
  int study, series;
  int x,y;

  int numFiles = 0, i=0;

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
    if ( fname.contains("thumbnails", Qt::CaseInsensitive) ||  fname.contains("dti", Qt::CaseInsensitive) ) 
  
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

      dcmDataSeries dsr(patient, protocol, study, series);
      // search to see if dsr is already in the list ...
      int idx = dseries.indexOf(dsr); 
      if ( idx == -1){
        dsr.files << fname;
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

        // preview ...
        const Uint16 *pixelData = NULL;
        unsigned char *prvu = new unsigned char[4*x*y];
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

  patientTreeWidget->clear();
  QList<QTreeWidgetItem*> _patientList;
  QList<QDomElement> domList;

  QDomElement root = patientTreeWidget->initNewRoot();
  QDomDocument* doc = patientTreeWidget->getDocument();
  // insert all patients as top level widgetitems ...
  for (i=0; i<patients.size(); i++) {
    QTreeWidgetItem *pat = new QTreeWidgetItem(patientTreeWidget);
    pat->setText(0, patients[i]);
    pat->setIcon(0, QIcon(QPixmap(":/res/images/patient.png")));
    // pat->setText(1, QString("dummy"));
    _patientList.append(pat);
	// DOM Stuff ...
	QDomElement tag = doc->createElement("patient");
	tag.setAttribute("folded", "yes");
	// tag.setAttribute("name", patients[i]);
    root.appendChild(tag);
	QDomElement tag2 = doc->createElement("name");
	tag.appendChild(tag2);
	QDomText t = doc->createTextNode(patients[i]);
    tag2.appendChild(t);
	domList.append(tag);
  }

  int frameCnt=1; 
  QDomElement lastFrame;
  QTreeWidgetItem *lastItem;

  QDir tnDir(directory);
  tnDir.cd("thumbnails");

  // now enter all the protocols ...
  for (i=0; i<dseries.size(); i++) {
	  // check if this in part of a multiframe dataset ...
	  if (i && ( dseries[i-1].patient == dseries[i].patient ) && ( dseries[i-1].protocol == dseries[i].protocol ) && ( dseries[i-1].x == dseries[i].x ) && ( dseries[i-1].y == dseries[i].y ) && ( dseries[i-1].files.size() == dseries[i].files.size() ) ) {
		  // simply append to existing dataset ...
		  frameCnt++;
		  // for tree View ...
		  lastItem->setText(5, QString::number(frameCnt));
		  // for DOM ...
		  lastFrame.setAttribute("count", frameCnt);
		  QFileInfo fi(dseries[i].files[0]);
		  QDomElement fr = doc->createElement("frame");
		  fr.setAttribute("directory", fi.dir().absolutePath() );
		  for (int f=0; f<dseries[i].files.size(); f++) {
			  QDomElement fn = doc->createElement("file");
			  QFileInfo fi2(dseries[i].files[f]);
			  fn.setAttribute("name", fi2.fileName());
			  fr.appendChild(fn);
		  }
		  lastFrame.appendChild(fr);
	  } else {
		  // reset frameCnt ...
		  frameCnt = 1;
		  int _pid = patients.indexOf(dseries[i].patient);
		  QTreeWidgetItem *srs = new QTreeWidgetItem( _patientList[_pid] );
		  srs->setText(0, dseries[i].patient);
		  srs->setText(1, dseries[i].protocol);
		  srs->setText(2, QString::number(dseries[i].study));
		  srs->setText(3, QString::number(dseries[i].series));
		  srs->setText(4, QString("%1x%2x%3").arg(dseries[i].x).arg(dseries[i].y).arg(dseries[i].files.size()));
		  srs->setText(5, QString::number(frameCnt));
		  // show a preview as the icon ...
		  QImage prv( dseries[i].preview, dseries[i].x, dseries[i].y, QImage::Format_RGB32 );
		  srs->setData(0, Qt::WhatsThisRole, QPixmap::fromImage(prv) );

                  // decide on filename ... and save it to the thumbnails
                  // folder ...
                  QStringList tnFilt;
                  tnFilt += "tn.*.png";
                  QStringList tnList = tnDir.entryList(tnFilt, QDir::Files, QDir::Name);
                  int tnCnt = tnList.size();  // @bug ... larger number could exist ...
                  QString tnName = QString().sprintf("tn.%06d.png", tnCnt);

                  // now save the thumbnail ....
                  prv.save( tnDir.absoluteFilePath(tnName), "PNG", 75 );

                  // free memory ...
                  delete [] dseries[i].preview;


		  srs->setIcon(0, QIcon(QPixmap(":/res/images/data.png"))); 
		  lastItem = srs;
		  // srs->setIcon(0, QIcon(QPixmap::fromImage(prv ))); 
		  // ui.m_lblPreview->setPixmap( QPixmap::fromImage( prv ) );

		  // now for the DOM entries ...
		  QDomElement tag = doc->createElement("dataset");
		  tag.setAttribute("type", "None");
		  // protocol ...
		  QDomElement proto = doc->createElement("protocol");
		  QDomText t = doc->createTextNode(dseries[i].protocol);
		  proto.appendChild(t);
		  tag.appendChild(proto);
		  // study ...
		  QDomElement study = doc->createElement("study");
		  QDomText t1 = doc->createTextNode( QString::number(dseries[i].study));
		  study.appendChild(t1);
		  tag.appendChild(study);
		  // study ...
		  QDomElement series = doc->createElement("series");
		  QDomText t2 = doc->createTextNode( QString::number(dseries[i].series));
		  series.appendChild(t2);
		  tag.appendChild(series);
		  // slices .. for now .. all dims later ...
		  QDomElement dims = doc->createElement("dimensions");
		  QDomText t3 = doc->createTextNode( QString("%1,%2,%3").arg(dseries[i].x).arg(dseries[i].y).arg(dseries[i].files.size()) );
		  dims.appendChild(t3);
		  tag.appendChild(dims);
		  // frames ...
		  QDomElement frames = doc->createElement("frames");
		  frames.setAttribute("count", frameCnt);
		  tag.appendChild(frames);
                  // thumbnail ...
                  QDomElement tnail = doc->createElement("thumbnail");
		  QDomText t4 = doc->createTextNode(tnName);
		  tnail.appendChild(t4);
                  frames.appendChild(tnail);

		  // for (int f=0; f<1; f++) {
		  QFileInfo fi(dseries[i].files[0]);
		  QDomElement fr = doc->createElement("frame");
		  fr.setAttribute("directory", fi.dir().absolutePath() );
		  for (int f=0; f<dseries[i].files.size(); f++) {
			  QDomElement fn = doc->createElement("file");
			  QFileInfo fi2(dseries[i].files[f]);
			  fn.setAttribute("name", fi2.fileName());
			  fr.appendChild(fn);
		  }
		  frames.appendChild(fr);
		  lastFrame = frames;
		  domList[_pid].appendChild(tag);
	  }
  }
  return true;
}
//
void patientBrowser::on_deleteButton_clicked()
{
}
//
void patientBrowser::on_importButton_clicked()
{
}
//
void patientBrowser::on_cdImportButton_clicked()
{
}
//
void patientBrowser::on_editPatientButton_clicked()
{
}
//
void patientBrowser::on_viewButton_clicked()
{
}
//
void patientBrowser::on_reloadButton_clicked() 
{
  rescanDataDirectory();
  // now write to a file ...
  QFile file(directory.absoluteFilePath("cardiac.xml"));
  file.open(QFile::WriteOnly | QFile::Text);
  patientTreeWidget->write(&file);
}
//
void patientBrowser::on_databaseButton_clicked()
{
}
//
void patientBrowser::on_exportButton_clicked()
{
}
//
void patientBrowser::on_anonymizeButton_clicked()
{
}
//
void patientBrowser::on_hdImportButton_clicked()
{
}
//

// Private utility functions ...
void patientBrowser::addFilesFromDir( const QString& dcmDir, QStringList& fileList) {
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

void patientBrowser::scanDir(const QString& dcmDir, QStringList& fileList) {
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

void patientBrowser::resizeCols() {
  for (int i=0; i<5; i++)
    patientTreeWidget->resizeColumnToContents(i);
}

void patientBrowser::updatePreview() {
  QList<QTreeWidgetItem *> items = patientTreeWidget->selectedItems();
  // only update on multiple items ...
  if (items.size() == 1) {
    QVariant pix = items[0]->data(0, Qt::WhatsThisRole);
    // update the preview ..
    previewLabel->setPixmap(pix.value<QPixmap>());
  }
}
