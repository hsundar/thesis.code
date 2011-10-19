#ifndef COLORSWATCH_H
#define COLORSWATCH_H

#include <QDockWidget>

class QTabWidget;
class QToolBox;
class QCheckBox;
class QSlider;
class QLabel;
class QRadioButton;
class QComboBox;

class tfEditor;
class viewer3d;

class ColorSwatch : public QDockWidget
{
  Q_OBJECT


  public:
    ColorSwatch(const QString &colorName, QWidget *parent = 0, Qt::WFlags flags = 0);
    void setViewer (viewer3d *view) { m_viewer = view; };

    // set labels ...
    void setVolDims(QString dims);
    void setVolSpacing(QString sp);
    void setPatientName(QString pn);

    void showTab(int i);
    void setMaxDepth(int i);

    protected slots:
      void sldrLowChanged(int);
    void sldrHighChanged(int);
    void sldrIsoChanged(int);
    void tfChanged();
    void isoApplyAllToggled(bool);
    void loadFA();
    void loadPD();
    void fiberTrack();

  protected:
    QCheckBox		*m_chkBBox, *m_chkClip, *m_chkVolume, *m_chkOctree;
    QCheckBox		*m_chkJac, *m_chkFibers, *m_chkIso, *m_chkDefFld;
    QComboBox		*m_FiberColorCombo;
    QSlider		*m_sldrOctHigh, *m_sldrOctLow, *m_sldrIso, *m_sldrFiberDensity;
    QLabel		*m_lblPatientName, *m_lblVolDims, *m_lblVolSpacing, *m_lblOctHigh, *m_lblOctLow, *m_lblIso;
    QTabWidget		*m_tabWidget;
    // QToolBox		*m_tabWidget;
    QRadioButton	*m_isoApplyAllBut;
  private:
    viewer3d	*m_viewer;
    tfEditor	*m_tfe;

    private slots:
};

#endif
