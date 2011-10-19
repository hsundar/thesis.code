#include "tfe.h"
#include "hoverpoints.h"
#include <iostream>
#include <fstream>

ShadeWidget::ShadeWidget(ShadeType type, QWidget *parent)
    : QWidget(parent), m_shade_type(type), m_alpha_gradient(QLinearGradient(0, 0, 0, 0))
{
    // Checkers background
    if (m_shade_type == ARGBShade) {
        QPixmap pm(20, 20);
        QPainter pmp(&pm);
        pmp.fillRect(0, 0, 10, 10, Qt::lightGray);
        pmp.fillRect(10, 10, 10, 10, Qt::lightGray);
        pmp.fillRect(0, 10, 10, 10, Qt::darkGray);
        pmp.fillRect(10, 0, 10, 10, Qt::darkGray);
        pmp.end();
        QPalette pal = palette();
        pal.setBrush(backgroundRole(), QBrush(pm));
        setAutoFillBackground(true);
        setPalette(pal);

    } else {
      setAttribute(Qt::WA_NoBackground);

    }

    QPolygonF points;
    points << QPointF(0.f, sizeHint().height())
      << QPointF(sizeHint().width(), 0.f);

    m_hoverPoints = new HoverPoints(this, HoverPoints::CircleShape);
    //     m_hoverPoints->setConnectionType(HoverPoints::LineConnection);
    m_hoverPoints->setPoints(points);
    m_hoverPoints->setPointLock(0, HoverPoints::LockToLeft);
    m_hoverPoints->setPointLock(1, HoverPoints::LockToRight);
    m_hoverPoints->setSortType(HoverPoints::XSort);

    // setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

    connect(m_hoverPoints, SIGNAL(pointsChanged(const QPolygonF &)), this, SIGNAL(colorsChanged()));
}


QPolygonF ShadeWidget::points() const
{
  return m_hoverPoints->points();
}


uint ShadeWidget::colorAt(int x)
{
  generateShade();

  QPolygonF pts = m_hoverPoints->points();
  for (int i=1; i < pts.size(); ++i) {
    if (pts.at(i-1).x() <= x && pts.at(i).x() >= x) {
      QLineF l(pts.at(i-1), pts.at(i));
      l.setLength(l.length() * ((x - l.x1()) / l.dx()));
      return m_shade.pixel(qRound(qMin(l.x2(), (qreal(m_shade.width() - 1)))),
          qRound(qMin(l.y2(), qreal(m_shade.height() - 1))));
    }
  }
  return 0;
}


void ShadeWidget::setGradientStops(const QGradientStops &stops)
{
  if (m_shade_type == ARGBShade) {
    m_alpha_gradient = QLinearGradient(0, 0, width(), 0);

    for (int i=0; i<stops.size(); ++i) {
      QColor c = stops.at(i).second;
      m_alpha_gradient.setColorAt(stops.at(i).first, QColor(c.red(), c.green(), c.blue()));
    }

    m_shade = QImage();
    generateShade();
    update();
  }
}

void ShadeWidget::paintEvent(QPaintEvent *)
{
  generateShade();

  QPainter p(this);
  p.drawImage(0, 0, m_shade);

  p.setPen(QColor(146, 146, 146));
  p.drawRect(0, 0, width() - 1, height() - 1);
}


void ShadeWidget::generateShade()
{
  if (m_shade.isNull() || m_shade.size() != size()) {

    if (m_shade_type == ARGBShade) {
      m_shade = QImage(size(), QImage::Format_ARGB32_Premultiplied);
      m_shade.fill(0);

      QPainter p(&m_shade);
      p.fillRect(rect(), m_alpha_gradient);

      p.setCompositionMode(QPainter::CompositionMode_DestinationIn);
      QLinearGradient fade(0, 0, 0, height());
      fade.setColorAt(0, QColor(0, 0, 0, 0));
      fade.setColorAt(1, QColor(0, 0, 0, 255));
      p.fillRect(rect(), fade);

    } else {
      m_shade = QImage(size(), QImage::Format_RGB32);
      QLinearGradient shade(0, 0, 0, height());
      shade.setColorAt(1, Qt::black);

      if (m_shade_type == RedShade)
        shade.setColorAt(0, Qt::red);
      else if (m_shade_type == GreenShade)
        shade.setColorAt(0, Qt::green);
      else
        shade.setColorAt(0, Qt::blue);

      QPainter p(&m_shade);
      p.fillRect(rect(), shade);
    }
  }
}


  tfEditor::tfEditor(QWidget *parent)
: QWidget(parent)
{
  // QGroupBox *editorGroup = new QGroupBox(this);
  // editorGroup->setAttribute(Qt::WA_ContentsPropagated);
  // editorGroup->setTitle("Transfer Function");
  // m_editor = new GradientEditor(editorGroup);

  QVBoxLayout *vbox = new QVBoxLayout(this);
  vbox->setSpacing(1);
  vbox->setMargin(1);

  m_red_shade = new ShadeWidget(ShadeWidget::RedShade, this);
  m_green_shade = new ShadeWidget(ShadeWidget::GreenShade, this);
  m_blue_shade = new ShadeWidget(ShadeWidget::BlueShade, this);
  m_alpha_shade = new ShadeWidget(ShadeWidget::ARGBShade, this);

  // m_red_shade->setMinimumWidth(200);
  // vbox->addWidget(editorGroup);
  vbox->addWidget(m_red_shade);
  vbox->addWidget(m_green_shade);
  vbox->addWidget(m_blue_shade);
  vbox->addWidget(m_alpha_shade);

  connect(m_red_shade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
  connect(m_green_shade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
  connect(m_blue_shade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
  connect(m_alpha_shade, SIGNAL(colorsChanged()), this, SLOT(pointsUpdated()));
}


inline static bool x_less_than(const QPointF &p1, const QPointF &p2)
{
  return p1.x() < p2.x();
}


void tfEditor::pointsUpdated()
{
  double w = m_alpha_shade->width();

  QGradientStops stops;

  QPolygonF points;

  points += m_red_shade->points();
  points += m_green_shade->points();
  points += m_blue_shade->points();
  points += m_alpha_shade->points();

  qSort(points.begin(), points.end(), x_less_than);

  for (int i=0; i<points.size(); ++i) {
    double x = int(points.at(i).x());
    if (i < points.size() - 1 && x == points.at(i+1).x())
      continue;
    QColor color((0x00ff0000 & m_red_shade->colorAt(int(x))) >> 16,
        (0x0000ff00 & m_green_shade->colorAt(int(x))) >> 8,
        (0x000000ff & m_blue_shade->colorAt(int(x))),
        (0xff000000 & m_alpha_shade->colorAt(int(x))) >> 24);

    if (x / w > 1)
      return;

    stops << QGradientStop(x / w, color);
  }

  m_alpha_shade->setGradientStops(stops);

  // emit gradientStopsChanged(stops);
  emit colormapChanged();
}

static void set_shade_points(const QPolygonF &points, ShadeWidget *shade)
{
  shade->hoverPoints()->setPoints(points);
  shade->hoverPoints()->setPointLock(0, HoverPoints::LockToLeft);
  shade->hoverPoints()->setPointLock(points.size() - 1, HoverPoints::LockToRight);
  shade->update();
}

void tfEditor::setGradientStops(const QGradientStops &stops)
{
  QPolygonF pts_red, pts_green, pts_blue, pts_alpha;

  double h_red = m_red_shade->height();
  double h_green = m_green_shade->height();
  double h_blue = m_blue_shade->height();
  double h_alpha = m_alpha_shade->height();

  for (int i=0; i<stops.size(); ++i) {
    double pos = stops.at(i).first;
    QRgb color = stops.at(i).second.rgba();
    pts_red << QPointF(pos * m_red_shade->width(), h_red - qRed(color) * h_red / 255);
    pts_green << QPointF(pos * m_green_shade->width(), h_green - qGreen(color) * h_green / 255);
    pts_blue << QPointF(pos * m_blue_shade->width(), h_blue - qBlue(color) * h_blue / 255);
    pts_alpha << QPointF(pos * m_alpha_shade->width(), h_alpha - qAlpha(color) * h_alpha / 255);
  }

  set_shade_points(pts_red, m_red_shade);
  set_shade_points(pts_green, m_green_shade);
  set_shade_points(pts_blue, m_blue_shade);
  set_shade_points(pts_alpha, m_alpha_shade);
  emit colormapChanged();
}

void tfEditor::updateColormap(unsigned char *cmap, int sz) {
  // get color from m_alpha_shade ... has rgba color ...
  QPolygonF points = m_alpha_shade->points();

  qSort(points.begin(), points.end(), x_less_than);

  float fac = ((float)sz)/points.at(points.size()-1).x();

  // std::ofstream out ("colors.txt");

  for (int i=0; i<points.size()-1; i++) {
    double x1 = int(points.at(i).x());
    double x2 = int(points.at(i+1).x());

    // if (i < points.size() - 1 && x1 == points.at(i+1).x())
    //	continue;

    QColor color1((0x00ff0000 & m_red_shade->colorAt(int(x1))) >> 16,
        (0x0000ff00 & m_green_shade->colorAt(int(x1))) >> 8,
        (0x000000ff & m_blue_shade->colorAt(int(x1))),
        (0xff000000 & m_alpha_shade->colorAt(int(x1))) >> 24);
    QColor color2((0x00ff0000 & m_red_shade->colorAt(int(x2))) >> 16,
        (0x0000ff00 & m_green_shade->colorAt(int(x2))) >> 8,
        (0x000000ff & m_blue_shade->colorAt(int(x2))),
        (0xff000000 & m_alpha_shade->colorAt(int(x2))) >> 24);

    // out << i << " " << color1.red() << " " << color1.green() << " " << color1.blue() << " " << color1.alpha() << std::endl;
    // out << i << " " << color2.red() << " " << color2.green() << " " << color2.blue() << " " << color2.alpha() << std::endl;
    int ind1 = (int)(x1*fac); int ind2 = (int)(x2*fac);
    int diff = ind2-ind1;
    for (int j=ind1; j<ind2; j++) {
      cmap[4*j] = color1.red() + (j-ind1)*(color2.red() - color1.red())/diff;
      cmap[4*j+1] = color1.green() + (j-ind1)*(color2.green() - color1.green())/diff;
      cmap[4*j+2] = color1.blue() + (j-ind1)*(color2.blue() - color1.blue())/diff;
      // Inventor makes alpha = 0.0 = opaque ...
      cmap[4*j+3] = color1.alpha() + (j-ind1)*(color2.alpha() - color1.alpha())/diff ;
    }
  }
  // out.close();
  // cmap[0] = 0; cmap[1]=0; cmap[2] = 0; cmap[3] = 255;
  // std::cout << "COLORMAP " << (int)cmap[0] << " " << (int)cmap[1] << " " << (int)cmap[2] << " " << (int)cmap[3] << std::endl;
}

void tfEditor::reset() {
  QGradientStops stops;

  stops << QGradientStop(0.00, QColor::fromRgba(0xff000000));
  stops << QGradientStop(1.00, QColor::fromRgba(0x00ffffff));

  setGradientStops(stops);
}

