#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDir>
#include <QStringList>

class QToolBar;
class QAction;
class QMenu;
class viewer3d;
class QDockWidget;
class toolDock;
class patientBrowser;

class MainWindow : public QMainWindow
{
  Q_OBJECT

    QToolBar *toolbar;
  // QMenu *cardiacMenu;

  public:
  MainWindow(QWidget *parent = 0, Qt::WFlags flags = 0);

  viewer3d* getViewer() {return m_viewer;};
  QDockWidget* getDock() {return m_dock;};
  toolDock* getToolDock() {return m_tools;};

  public slots:
    void openRecentFile();
  void actionTriggered(QAction *action);
  void loadVolume();
  void loadField();
  void unloadVolume();
  void loadOctree();
  void generateOctree();
  void showDcmDlg();
  void showSimDlg();
  void showPatientBrowser();
  void about();
  void aboutPlugins();
  void changeBackground();
  void alignVolumes();
	void animateLandmarks();

  protected:
  void closeEvent(QCloseEvent *event);

  // variables ...
  viewer3d	*m_viewer;
  QDockWidget	*m_dock;
  toolDock    *m_tools;

  QMenu *recentFilesMenu;
  enum { MaxRecentFiles = 5 };
  QAction *recentFileActs[MaxRecentFiles];
  QAction *loadVolAct;
  QAction *loadDefAct;
  QAction *unloadVolAct;
  QAction *loadOctAct;
  QAction *genOctAct;
  QAction *alignAct;
  QAction *changeBGAct;
  QAction *dicomAct;
  QAction *patBrowseAct;
  QAction *aboutPluginsAct;
  QAction *separatorAct;
  QAction *simAct;
  QAction *aboutAct;

	QAction *animateLandmarksAct;
private:
  void setupToolBar();
  void setupActions();
  void setupMenuBar();
  void setupDockWidgets();
  void loadPlugins();

	void readSettings();
  void writeSettings();

  void setCurrentFile(const QString &fileName);
  void updateRecentFileActions();

  QDir			pluginsDir;
  QDir			databaseDir;
  QStringList		pluginFileNames;

  patientBrowser	*pBrowser;
};

#endif
