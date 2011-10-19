#ifndef __VIEWER_3D_H_
#define __VIEWER_3D_H_

#include <QtCore>
#include <QtGui>

class Volume;

class SoSeparator;
class SoVolumeData; 	
class SoTransferFunction;
class SoVolumeRender;
class BoundingBox;		
class SoClipPlaneManip;
class SoSwitch;
class SoBlinker;
class SbVec2s;
class MarchingCubes;
class SoMaterial;
class FiberTracker;
class SoVectorFieldViz;
class MyExaminerViewer;

#include "Point.h"
#include "Field.h"

class viewer3d : public QWidget {
  Q_OBJECT
  public:
    enum playType {PLAY, PAUSE, REW, FWD};
    viewer3d( QWidget * parent = 0, Qt::WFlags f = 0 );

    void loadVolume (const QString &filename);
    void loadOctree (const QString &filename);
    void loadField (const QString &filename);
    void loadPoints(const QString &filename);

		void animateLandmarks(const QString &dirname);

    void unloadVolume();

    void alignVolumes();

    void setOctThreshold(int, int);

    unsigned char* getColormap() { return m_ucColors; };
    void colormapChanged();
    void redrawDeformationField();

    void generateOctree();
    void controlSequencer(playType mode);
    void setSequencerSpeed(double speed);
    void changeIsoLevel(int lev);
    void updateIsoAll(int lev);

    // DTI / Fiber track...
    void loadFA(const QString &fname);
    void loadPD(const QString &fname);
    void fiberTrack(int density, int color);

    // set/get BG color
    QColor getBackgroundColor();
    void setBackgroundColor(QColor col);

    void toggleClipPlane();
    void toggleBoundingBox();
    void toggleVolume();
    void toggleOctree();

    QList<Volume*> getVolumes();
  protected slots:
    void playSeq();
    void pauseSeq();
    void rewSeq();
    void fwdSeq();

  protected:
    SoSeparator* addIsosurface(Volume*);

    // Variables ...
    MyExaminerViewer 	*m_eView;

    // Data
    QList<Volume*>		m_volume;
    //Point			m_volSize;
    //Point			m_volSpacing;

    // Inventor ...
    SoSeparator		*m_sceneGraph;
    SoSeparator		*m_volGroup;
    SoSeparator		*m_geoGroup;

    // Voleon Stuff ...
    QList<SoVolumeData*> 	m_soVolumeData;
    SoTransferFunction      *m_soTF;
    SoVolumeRender 			*m_soVolRen;

    // Isosurface ...
    QList<MarchingCubes*>   m_marchers;
    QList<SoMaterial*>      m_isoMaterial;

    // deformation field visualization ...
    SoVectorFieldViz		*m_defViz;
    QList< Field<float, 3>* >			m_defField;

    // DTI, FIber tracking ...
    FiberTracker			*m_tracker;
    Field<float, 3>			*m_pPD;
    Field<float, 1>			*m_pFA;

    // Decorator stuff ...
    BoundingBox		*m_soBoundingBox;
    SoClipPlaneManip	*m_soClipPlane;

    // Switches ...
    SoSwitch		*m_soClipSwitch;
    SoSwitch		*m_soBBoxSwitch;
    SoBlinker		*m_soVolSwitch; // changed switch to blinker ...
    SoSwitch		*m_soOctSwitch;

    // Octree
    int 			m_iDepth;

    // flags
    bool 			m_bHaveVolume;
    bool 			m_bHaveOctree;
		bool 			m_bHaveField;

    // colormap
    unsigned char		*m_ucColors;

		// Landmarks
		QList< QList <Point>> m_ptLandmarks;
};

#endif
