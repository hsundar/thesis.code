#include "mainwindow.h"
#include "colorswatch.h"
#include "tooldock.h"
#include "plugindialog.h"

#include <QAction>
#include <QLayout>
#include <QMenu>
#include <QMenuBar>
#include <QStatusBar>
#include <QTextEdit>
#include <QToolBar>

#include "viewer3d.h"
#include "dcmDialog.h"
#include "patientBrowser.h"

#include <fstream>
#include <iostream>

MainWindow::MainWindow(QWidget *parent, Qt::WFlags flags)
    : QMainWindow(parent, flags)
{
    pBrowser = NULL;

    setObjectName("MainWindow");
    setWindowTitle("Cardiac Viewer");

    // qDebug("Creating new viewer3d");
    //
    m_viewer = new viewer3d(this);
    m_viewer->setMinimumSize(600, 600);

    setupDockWidgets();
    setupActions();
    setupToolBar();
    setupMenuBar();
    loadPlugins();

    readSettings();

    setCentralWidget(m_viewer);
    statusBar()->showMessage(tr("Status Bar"));
}

void MainWindow::actionTriggered(QAction *action)
{
	qDebug("action '%s' triggered", action->text().toLocal8Bit().data());
}

void MainWindow::setupToolBar()
{
  toolbar = addToolBar(tr("&Tools"));
  toolbar->setAllowedAreas(Qt::TopToolBarArea | Qt::LeftToolBarArea);
  toolbar->addAction(loadVolAct);
  toolbar->addAction(loadOctAct);
  toolbar->addAction(dicomAct);
  toolbar->addAction(patBrowseAct);
  toolbar->addAction(changeBGAct);

}

void MainWindow::setupActions() {

    for (int i = 0; i < MaxRecentFiles; ++i) {
        recentFileActs[i] = new QAction(this);
        recentFileActs[i]->setVisible(false);
        connect(recentFileActs[i], SIGNAL(triggered()),
                this, SLOT(openRecentFile()));
    }

    loadVolAct = new QAction(tr("Load &Volume"), this);
    loadVolAct->setStatusTip(tr("Load a volume into the viewer"));
    loadVolAct->setShortcut(QKeySequence(tr("Ctrl+V")));
    loadVolAct->setIcon(QIcon(QPixmap(":/res/images/loadvol.png")));
    connect(loadVolAct, SIGNAL(triggered()), this, SLOT(loadVolume()) );

    loadOctAct = new QAction(tr("Load &Octree"), this);
    loadOctAct->setStatusTip(tr("Load an Octree into the viewer"));
    loadOctAct->setShortcut(QKeySequence(tr("Ctrl+O")));
    loadOctAct->setIcon(QIcon(QPixmap(":/res/images/loadoct.png")));
    connect(loadOctAct, SIGNAL(triggered()), this, SLOT(loadOctree())) ;

    alignAct = new QAction(tr("&Align Volumes"), this);
    alignAct->setStatusTip(tr("Align two volumes by performing Rigid Registration"));
    alignAct->setShortcut(QKeySequence(tr("Ctrl+A")));
    alignAct->setIcon(QIcon(QPixmap(":/res/images/registration.png")));
    connect(alignAct, SIGNAL(triggered()), this, SLOT(alignVolumes()));

		animateLandmarksAct = new QAction(tr("Animate Landmarks"), this);
		animateLandmarksAct->setStatusTip(tr("Animate landmarks using the deformation field"));
		connect(animateLandmarksAct, SIGNAL(triggered()), SLOT(animateLandmarks())) ;


	loadDefAct = new QAction(tr("Load Vector &Field"), this);
	loadDefAct->setStatusTip(tr("Load a deformation field into the viewer"));
	loadDefAct->setShortcut(QKeySequence(tr("Ctrl+F")));
	loadDefAct->setIcon(QIcon(QPixmap(":/res/images/loadoct.png")));
    connect(loadDefAct, SIGNAL(triggered()), this, SLOT(loadField())) ;

    genOctAct = new QAction(tr("&Generate Octree"), this);
    genOctAct->setStatusTip(tr("Generate octree based on the loaded volume dataset"));
    genOctAct->setShortcut(QKeySequence(tr("Ctrl+G")));
    genOctAct->setIcon(QIcon(QPixmap(":/res/images/octree.png")));
    connect(genOctAct, SIGNAL(triggered()), this, SLOT(generateOctree()));
    
    changeBGAct = new QAction(tr("Change &Background Color"), this);
    changeBGAct->setStatusTip(tr("Change the background color for the Viewer window"));
    changeBGAct->setShortcut(QKeySequence(tr("Ctrl+C")));
    changeBGAct->setIcon(QIcon(QPixmap(":/res/images/display.png")));
    connect(changeBGAct, SIGNAL(triggered()), this, SLOT(changeBackground()));

	simAct = new QAction(tr("Similarity Measures"), this);
	simAct->setStatusTip(tr("Open up the dialog to evaluate similarity measures and perform affine registration"));
	simAct->setIcon(QIcon(QPixmap(":/res/images/registration.png")));
	connect(simAct, SIGNAL(triggered()), this, SLOT(showSimDlg()));

    dicomAct = new QAction(tr("Open Dicom"), this);
    dicomAct->setStatusTip(tr("Open the Dicom import dialog"));
    dicomAct->setShortcut(QKeySequence(tr("Ctrl+D")));
    dicomAct->setIcon(QIcon(QPixmap(":/res/images/dicom.png")));
    connect(dicomAct, SIGNAL(triggered()), this, SLOT(showDcmDlg()));

		patBrowseAct = new QAction(tr("Open &Patient Browser"), this);
		patBrowseAct->setStatusTip(tr("Open the patient browse"));
		patBrowseAct->setShortcut(QKeySequence(tr("Ctrl+P")));
		patBrowseAct->setIcon(QIcon(QPixmap(":/res/images/patient.png")));
		connect(patBrowseAct, SIGNAL(triggered()), this, SLOT(showPatientBrowser())) ;

    aboutPluginsAct = new QAction(tr("About &Plugins"), this);
    connect(aboutPluginsAct, SIGNAL(triggered()), this, SLOT(aboutPlugins()));

    aboutAct = new QAction(tr("&About"), this);
    aboutAct->setStatusTip(tr("About Cardiac Viewer"));
    connect(aboutAct, SIGNAL(triggered()), this, SLOT(about()));
}
void MainWindow::setupMenuBar()
{
  QMenu *menu = menuBar()->addMenu(tr("&File"));
  menu->addAction(loadVolAct);
  menu->addAction(loadDefAct);
  menu->addAction(loadOctAct);
  separatorAct = menu->addSeparator();
  for (int i = 0; i < MaxRecentFiles; ++i)
    menu->addAction(recentFileActs[i]);
  menu->addSeparator();
  menu->addAction(tr("&Quit"), this, SLOT(close()));
  updateRecentFileActions();

  menuBar()->addSeparator();

  QMenu *utils = menuBar()->addMenu(tr("&Utils"));
  utils->addAction(genOctAct);
  utils->addAction(dicomAct);
  utils->addAction(patBrowseAct);
  utils->addAction(changeBGAct);
  // utils->addAction(simAct);
  utils->addAction(alignAct);
	utils->addAction(animateLandmarksAct);

  QMenu *view = menuBar()->addMenu(tr("&View"));
  view->addAction(m_dock->toggleViewAction()); 

  QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutPluginsAct);
  helpMenu->addAction(aboutAct);
}

void MainWindow::setupDockWidgets()
{
  m_dock = new QDockWidget(tr("Toolchest"), this);
  m_dock->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea );
  m_dock->setMaximumWidth(300);
  m_tools = new toolDock( tr("Toolchest"), m_dock);
  m_dock->setWidget(m_tools);
  m_tools->setViewer(m_viewer);
  
  addDockWidget(Qt::RightDockWidgetArea, m_dock);
}

void MainWindow::openRecentFile()
{
  QAction *action = qobject_cast<QAction *>(sender());
  if (action) {
    QString fname = action->data().toString();
    if ( fname.endsWith("oct") )
      m_viewer->loadOctree(fname);
    if ( fname.endsWith("mhd") || fname.endsWith("hdr") )
      m_viewer->loadVolume(fname);

    setCurrentFile(fname);
    statusBar()->showMessage(tr("File loaded"), 2000);
  }
}

void MainWindow::loadVolume()
{
  QStringList files = QFileDialog::getOpenFileNames(
			this,
			"Choose a file",
			"/home",
			"MetaIO/Analyze Files (*.mhd *.mha *.hdr)");
 
  QStringListIterator i(files);
 
  // set up progress dialog for this ...
  int cnt=0;
  QProgressDialog progress("Loading Volumes...", "Abort Load", 0, files.size(), this);
  progress.setValue(cnt);

  while (i.hasNext()) {
    progress.setValue(cnt++);
    qApp->processEvents();
    QString fn = i.next();
    // std::cout << "Loading " << qPrintable(fn) << std::endl;
    m_viewer->loadVolume( fn );
    setCurrentFile(fn);
    statusBar()->showMessage(fn + " loaded", 2000);

    if (progress.wasCanceled())
      break;
    }
  progress.setValue(files.size());
  // if (files.size() > 3)
  //   m_viewer->controlSequencer(viewer3d::PLAY);
}

void MainWindow::loadOctree()
{
    QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "/home", "Octree Files (*.oct)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadOctree( fn );
	setCurrentFile(fn);
	statusBar()->showMessage(tr("Octree File loaded"), 2000);
    }
}

void MainWindow::loadField()
{
    /*QString fn = QFileDialog::getOpenFileName( this, "Choose a file", "/home", "MetaIO Files (*.mhd)" );
    if ( !fn.isEmpty() ) {
        m_viewer->loadField( fn );
	setCurrentFile(fn);
	statusBar()->showMessage(tr("Vector FIeld loaded"), 2000);
    }*/

		QStringList files = QFileDialog::getOpenFileNames(
			this,
			"Choose a file",
			"c:/home/data",
			"MetaIO/Analyze Files (*.mhd *.mha *.hdr)");

		QStringListIterator i(files);

		// set up progress dialog for this ...
		int cnt=0;
		QProgressDialog progress("Loading Fields...", "Abort Load", 0, files.size(), this);
		progress.setValue(cnt);

		while (i.hasNext()) {
			progress.setValue(cnt++);
			qApp->processEvents();
			QString fn = i.next();
			// std::cout << "Loading " << qPrintable(fn) << std::endl;
			m_viewer->loadField( fn );
			// setCurrentFile(fn);
			statusBar()->showMessage(fn + " loaded", 2000);

			if (progress.wasCanceled())
				break;
		}
		progress.setValue(files.size());
}

void MainWindow::setCurrentFile(const QString &fileName)
{
    QSettings settings("Hari Sundar", "Cardiac Viewer");
    QStringList files = settings.value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > MaxRecentFiles)
        files.removeLast();

    settings.setValue("recentFileList", files);

    foreach (QWidget *widget, QApplication::topLevelWidgets()) {
        MainWindow *mainWin = qobject_cast<MainWindow *>(widget);
        if (mainWin)
            mainWin->updateRecentFileActions();
    }
}

void MainWindow::updateRecentFileActions()
{
    QSettings settings("Hari Sundar", "Cardiac Viewer");
    QStringList files = settings.value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

    for (int i = 0; i < numRecentFiles; ++i) {
        QString text = tr("&%1 %2").arg(i + 1).arg(QFileInfo(files[i]).fileName());
        recentFileActs[i]->setText(text);
        recentFileActs[i]->setData(files[i]);
        recentFileActs[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
        recentFileActs[j]->setVisible(false);

    separatorAct->setVisible(numRecentFiles > 0);
}


void MainWindow::about()
{
   QMessageBox::about(this, tr("About Cardiac Viewer"),
            tr("Hari Sundar\nSection for Biomedical Image Analysis\nUniversity of Pennsylvania\n\nContact: hsundar@seas.upenn.edu\n\n"
	       "DISCLAIMER\n\n"   
	       "This Software is a research prototype and is not intended to be used in a clinical setting.The\n"
	       "authors will not be held responsible for errors or technical difficulties with the  functionality\n"
	       "of the software. By using this software,  you assume any and all responsibility and risk of use of\n"
	       "the software. The software is provided on an as is basis without warranties of any kind either\n"
	       "expressed or implied. The authors may make improvements to the software at any time with or\n" 
	       "without notice and will not guarantee the correction of software for errors reported.\n"
	       ));
}

void MainWindow::aboutPlugins()
{
    PluginDialog dialog(pluginsDir.path(), pluginFileNames, this);
    dialog.exec();
}

void MainWindow::loadPlugins()
{
    pluginsDir = QDir(qApp->applicationDirPath());

#if defined(Q_OS_WIN)
    if (pluginsDir.dirName().toLower() == "debug" || pluginsDir.dirName().toLower() == "release")
        pluginsDir.cdUp();
#elif defined(Q_OS_MAC)
    if (pluginsDir.dirName() == "MacOS") {
        pluginsDir.cdUp();
        pluginsDir.cdUp();
        pluginsDir.cdUp();
    }
#endif
    pluginsDir.cd("plugins");

    foreach (QString fileName, pluginsDir.entryList(QDir::Files)) {
        QPluginLoader loader(pluginsDir.absoluteFilePath(fileName));
        QObject *plugin = loader.instance();
        if (plugin) {
           /* RegInterface *iReg = qobject_cast<RegInterface *>(plugin);
            if (iReg)
                addToMenu(plugin, iReg->algos(), algoMenu,
                          SLOT(changeBrush()), brushActionGroup );

            ShapeInterface *iShape = qobject_cast<ShapeInterface *>(plugin);
            if (iShape)
                addToMenu(plugin, iShape->shapes(), shapesMenu,
                          SLOT(insertShape()));

            FilterInterface *iFilter = qobject_cast<FilterInterface *>(plugin);
            if (iFilter)
                addToMenu(plugin, iFilter->filters(), filterMenu,
                          SLOT(applyFilter()));
            */
            pluginFileNames += fileName;
        }
    }

    // brushMenu->setEnabled(!brushActionGroup->actions().isEmpty());
    // shapesMenu->setEnabled(!shapesMenu->actions().isEmpty());
    // filterMenu->setEnabled(!filterMenu->actions().isEmpty());
}

void MainWindow::generateOctree() {
  m_viewer->generateOctree();
}

void MainWindow::showDcmDlg() {
  dcmDialog dlg(this);
  dlg.exec();
}

void MainWindow::showSimDlg() {
}

void MainWindow::alignVolumes() {
  m_viewer->alignVolumes();
}

void MainWindow::showPatientBrowser() {
  if (!pBrowser) {
    pBrowser = new patientBrowser(this);
    pBrowser->initDataDirectory(databaseDir);
  }
  pBrowser->show();
  pBrowser->activateWindow();
}

void MainWindow::changeBackground() {
  // get current color from viewer ...
  QColor oldCol = m_viewer->getBackgroundColor();
  // pop up a color selection dialog ...
  QColor col = QColorDialog::getColor ( oldCol, this );

  // now set it
  m_viewer->setBackgroundColor(col);
}

void MainWindow::unloadVolume() {
  m_viewer->unloadVolume();
}

void MainWindow::readSettings() {
  QSettings settings("Hari Sundar", "Cardiac Viewer");

  QRect rect = settings.value("geometry", QRect(100,100, 800, 800)).toRect();

  move(rect.topLeft());
  resize(rect.size());

  QString str = settings.value("Directories/plugins").toString();

  if ( str == "" ) {
    QString s = QFileDialog::getExistingDirectory(
                    this,
                    "Please select the Plugins directory",
                    QDir::homePath(),
                    QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    str = s;
  }

  // std::cout << "Plugins dir: " << qPrintable(str) << std::endl;
  
  pluginsDir.cd(str);
  str = settings.value("Directories/database").toString();
  if ( str == "" ) {
    QString s = QFileDialog::getExistingDirectory(
                    this,
                    "Please select the Patient Database directory",
                    QDir::homePath(),
                    QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    str = s;
  }

  // std::cout << "Database dir: " << qPrintable(str) << std::endl;
  databaseDir.cd(str); 

  // check if directory exists for thumbnails .. else create it ...
  if (! databaseDir.exists("thumbnails") ) {
    // create the thumbnails dir ...
    databaseDir.mkdir("thumbnails");
  }

}

void MainWindow::writeSettings() {
  QSettings settings("Hari Sundar", "Cardiac Viewer");

  settings.setValue("geometry", geometry());
  settings.setValue("Directories/database", databaseDir.absolutePath());
  settings.setValue("Directories/plugins", pluginsDir.absolutePath());
}

void MainWindow::closeEvent(QCloseEvent *event) {
  writeSettings();
  event->accept();
}

void MainWindow::animateLandmarks() {
	// get directory to animate from
	QString s = QFileDialog::getExistingDirectory(
		this,
		"Please select the data directory directory",
		QDir::homePath(),
		QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

	m_viewer->animateLandmarks(s);
}