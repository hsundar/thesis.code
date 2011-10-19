#ifndef __TFE_H_
#define __TFE_H_

#include "arthurwidgets.h"

#include <QtGui>

class HoverPoints;

class ShadeWidget : public QWidget
{
    Q_OBJECT
public:
    enum ShadeType {
        RedShade,
        GreenShade,
        BlueShade,
        ARGBShade
    };

    ShadeWidget(ShadeType type, QWidget *parent);

    void setGradientStops(const QGradientStops &stops);

    void paintEvent(QPaintEvent *e);

    QSize sizeHint() const { return QSize(150, 40); }
    QPolygonF points() const;

    HoverPoints *hoverPoints() const { return m_hoverPoints; }

    uint colorAt(int x);

signals:
    void colorsChanged();

private:
    void generateShade();

    ShadeType m_shade_type;
    QImage m_shade;
    HoverPoints *m_hoverPoints;
    QLinearGradient m_alpha_gradient;
};

class tfEditor : public QWidget
{
    Q_OBJECT
public:
    tfEditor(QWidget *parent);

    void setGradientStops(const QGradientStops &stops);
    void updateColormap(unsigned char *cmap, int sz);

    void reset();
    
public slots:
    void pointsUpdated();

signals:
    // void gradientStopsChanged(const QGradientStops &stops);
    void colormapChanged();

private:
    ShadeWidget *m_red_shade;
    ShadeWidget *m_green_shade;
    ShadeWidget *m_blue_shade;
    ShadeWidget *m_alpha_shade;
};

#endif // __TFE_H_
