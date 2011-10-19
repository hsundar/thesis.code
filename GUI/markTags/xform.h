#ifndef XFORM_H
#define XFORM_H

#include "arthurwidgets.h"

#include <QPolygonF>
#include <QList>

class HoverPoints;
class QLineEdit;
class QSlider;
class QSpinBox;
class Volume;

#define LM_SIZE 6

class XFormView : public ArthurFrame
{
    Q_OBJECT

    Q_PROPERTY(double rotation READ rotation WRITE setRotation)
    Q_PROPERTY(double scale READ scale WRITE setScale)

public:
    XFormView(QWidget *parent);
    void paint(QPainter *);
    void drawVectorType(QPainter *painter);
    void drawGrid(QPainter *painter);
    void drawPixmapType(QPainter *painter);
    void drawTextType(QPainter *painter);
    QSize sizeHint() const { return QSize(800, 800); }

    void mousePressEvent(QMouseEvent *e);
    void resizeEvent(QResizeEvent *e);
    HoverPoints *hoverPoints() { return ctrl; }

    double scale() const { return m_scale; }
    double rotation() const { return m_rotation; }
    void setScale(double s);
    void setRotation(double r);

public slots:
    void updateCtrlPoints(const QPolygonF &);
    void updateGridPoints(const QPolygonF &);
    void changeRotation(int rotation);
    void changeScale(int scale);
    void changeSlice(int s);
    void changeSpacing(int s);

    void copyPrev();
    void copyNext();

    void setVectorType();
    void setPixmapType();
    void setTextType();
    void reset();

    void load();
    void loadLM();
    void saveLM();
      
    void enableHover(int flag);
		void enableGrid(int flag);
    void enablePath(int flag);
      
signals:
    void rotationChanged(int rotation);
    void scaleChanged(int scale);
    void slicesChanged(int slices);
    void enablePts(bool);

protected:
    void wheelEvent(QWheelEvent *);

    QPolygonF screenToImage(QPolygonF pts);
    QPolygonF imageToScreen(QPolygonF pts);
		
		inline QRectF pointBoundingRect(int i) const;

private:
    enum XFormType { VectorType, PixmapType, TextType };

    QPolygonF ctrlPoints;
    QPolygonF gridPoints;
    QList<QPolygonF> imagePoints;
    HoverPoints *ctrl;
    HoverPoints *grid;

    int gridSpacing;

    unsigned int m_uiX;
    unsigned int m_uiY;

    double m_rotation;
    double m_scale;
    XFormType type;
    QList<QPixmap> pixmap;

    int curr;

    int gridOpacity;

    Volume * m_volume;
		
    bool m_bHaveVolume;
		bool m_bDrawGrid;
    bool m_bDrawPath;
		
		int m_iLastCopyWasFromNext;
};

inline QRectF XFormView::pointBoundingRect(int i) const
{
	QPointF p = gridPoints.at(i);
	double w = LM_SIZE;
	double h = LM_SIZE;
	double x = p.x() - w / 2;
	double y = p.y() - h / 2;
	return QRectF(x, y, w, h);
}


class XFormWidget : public QWidget
{
    Q_OBJECT
public:
    XFormWidget(QWidget *parent);

public slots:
  void dcm();
  void setNumSlices(int n);
private:
    XFormView *view;

    // GUI items 
    QSlider   *sliceSlider;
    QSpinBox  *sliceSpin;
};

#endif // XFORM_H
