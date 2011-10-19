#include "xform.h"
#include "hoverpoints.h"

#include <QLayout>
#include <QPainter>
#include <QPainterPath>

#include "Volume.h"
#include "dcmDialog.h"

#include <iostream>
#include <fstream>

const int alpha = 55;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


XFormView::XFormView(QWidget *parent)
: ArthurFrame(parent) {
  setAttribute(Qt::WA_MouseTracking);
  // type = VectorType;
  type = TextType;
  m_rotation = 0.0;
  m_scale = 1.0;
  m_bHaveVolume = false;
  m_bDrawGrid = false;
  m_bDrawPath = false;
  gridSpacing = 4;
  gridOpacity = 100;
	m_iLastCopyWasFromNext = 0;

  // Control Points
  ctrl = new HoverPoints(this, HoverPoints::RectangleShape);
  ctrl->setConnectionType(HoverPoints::LineConnection);
  ctrl->setEditable(false);
  ctrl->setPointSize(QSize(10, 10));
  ctrl->setShapeBrush(QBrush(QColor(151, 0, 0, alpha)));
  ctrl->setShapePen(QPen(QColor(255, 100, 50, alpha)));
  ctrl->setConnectionPen(QPen(QColor(151, 0, 0, alpha), 0, Qt::DotLine, Qt::FlatCap, Qt::BevelJoin));
  ctrl->setBoundingRect(QRectF(0, 0, 800, 800));
  ctrlPoints << QPointF(400, 400) << QPointF(550, 400);
  ctrl->setPoints(ctrlPoints);
  connect(ctrl, SIGNAL(pointsChanged(const QPolygonF&)),
          this, SLOT(updateCtrlPoints(const QPolygonF &)));

  // Grid ...
  grid = new HoverPoints(this, HoverPoints::CircleShape);
  grid->setConnectionType(HoverPoints::NoConnection);
  grid->setEditable(false);
	grid->setResizable(false);
  grid->setPointSize(QSize(LM_SIZE, LM_SIZE));
  grid->setShapeBrush(QBrush(QColor(0, 151, 0, 155)));
  grid->setShapePen(QPen(QColor(100, 255, 50, 155)));
  grid->setConnectionPen(QPen(QColor(0, 151, 0, 155), 0, Qt::DotLine, Qt::FlatCap, Qt::BevelJoin));
  grid->setBoundingRect(QRectF(0, 0, 800, 800));

  connect(grid, SIGNAL(pointsChanged(const QPolygonF&)),
          this, SLOT(updateGridPoints(const QPolygonF &)));

  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  curr=0;

  connect (this, SIGNAL(enablePts(bool)), ctrl, SLOT(setEnabled(bool)));
}

void XFormView::enableHover(int fl) {
  emit enablePts(fl);
}

void XFormView::enableGrid(int fl) {
  m_bDrawGrid = fl;
  update();
}

void XFormView::enablePath(int fl) {
  m_bDrawPath = fl;
  update();
}

void XFormView::mousePressEvent(QMouseEvent* ev) {
  setDescriptionEnabled(false);

  // handle clicks ...

  if ((ev->button() == Qt::LeftButton) && ((ev->modifiers() & Qt::ShiftModifier) == Qt::ShiftModifier)) {
    // pick point ...
    // std::cout << "Picked point is " << ev->x() << " " << ev->y() << std::endl;
    gridPoints << ev->pos(); // QPointF(x(), y());
    imagePoints[curr] = screenToImage(gridPoints);
    grid->setPoints(gridPoints);
  }

  QPointF clickPos = ev->pos();
  if ( ev->button() == Qt::RightButton) {
    // find the point ...
    int index = -1;
    for (int i=0; i<gridPoints.size(); ++i) {
      QPainterPath path;
      path.addEllipse(pointBoundingRect(i));

      if (path.contains(clickPos)) {
        index = i;
        break;
      }
    }   
    // delete the point ...
    if (index >=0) {
      gridPoints.remove(index);
      imagePoints[curr] = screenToImage(gridPoints);
      grid->setPoints(gridPoints);
    }
  }

  update();
}

void XFormView::resizeEvent(QResizeEvent *ev) {
  ctrl->setBoundingRect(rect());
  grid->setBoundingRect(rect());
	
	// now reset the grid Points ...
	if (m_bHaveVolume) {
    gridPoints = imageToScreen(imagePoints[curr]);
    grid->setPoints(gridPoints);
  }
}

void XFormView::paint(QPainter *p) {
  p->save();
  p->setRenderHint(QPainter::Antialiasing);
  p->setRenderHint(QPainter::SmoothPixmapTransform);
  switch (type) {
  case VectorType:
    drawVectorType(p);
    break;
  case PixmapType:
    drawPixmapType(p);
    // drawGrid(p);
    break;
  case TextType:
    drawTextType(p);
    break;
  }
  p->restore();
}

void XFormView::updateCtrlPoints(const QPolygonF &points) {
  QPointF trans = points.at(0) - ctrlPoints.at(0);

  if (qAbs(points.at(0).x() - points.at(1).x()) < 10
      && qAbs(points.at(0).y() - points.at(1).y()) < 10)
    ctrl->setPoints(ctrlPoints);
  if (!trans.isNull()) {
    ctrlPoints[0] = points.at(0);
    ctrlPoints[1] += trans;
    ctrl->setPoints(ctrlPoints);

    if (m_bHaveVolume) {
      gridPoints = imageToScreen(imagePoints[curr]);
      grid->setPoints(gridPoints);
    }
  }
  ctrlPoints = points;

  QLineF line(ctrlPoints.at(0), ctrlPoints.at(1));
  m_rotation = line.angle(QLineF(0, 0, 1, 0));
  if (line.dy() < 0)
    m_rotation = 360 - m_rotation;

  if (trans.isNull())
    emit rotationChanged(int(m_rotation*10));
}


void XFormView::updateGridPoints(const QPolygonF &points) {
  // Points are in screen space ...
  // Should store them internally in image space ...
  gridPoints = points;
  if (m_bHaveVolume)
    imagePoints[curr] = screenToImage(gridPoints);
  grid->setPoints(gridPoints);
}


void XFormView::setVectorType() {
  type = VectorType;
  update();
}

void XFormView::setPixmapType() {
  type = PixmapType;
  update();
}

void XFormView::setTextType() {
  type = TextType;
  update();
}

void XFormView::changeRotation(int r) {
  setRotation(double(r)/10.0);
}

void XFormView::changeScale(int s) {
  setScale(double(s)/1000.0);
}

void XFormView::changeSpacing(int s) {
  gridSpacing = s;
  update();
}

void XFormView::changeSlice(int s) {
  curr = s;
  if (curr >= pixmap.size())
    curr = 0;
  // change the gridPoints ...
  if (m_bHaveVolume) {
    gridPoints = imageToScreen(imagePoints[curr]);
    grid->setPoints(gridPoints);
  }
  update();
}

void XFormView::setScale(double s) {
  m_scale = s;
  if (m_bHaveVolume) {
    gridPoints = imageToScreen(imagePoints[curr]);
    grid->setPoints(gridPoints);
  }
  update();
}

void XFormView::setRotation(double r) {
  double old_rot = m_rotation;
  m_rotation = r;

  QPointF center(ctrl->points().at(0));
  QMatrix m;
  m.translate(center.x(), center.y());
  m.rotate(m_rotation - old_rot);
  m.translate(-center.x(), -center.y());
  ctrl->setPoints(ctrl->points() * m);

  if (m_bHaveVolume) {
    gridPoints = imageToScreen(imagePoints[curr]);
    grid->setPoints(gridPoints);
  }

  update();
}

void XFormView::wheelEvent(QWheelEvent *e) {
  m_scale += e->delta()/600.0;
  m_scale = qMax(.1, qMin(12.0, m_scale));
  emit scaleChanged(int(m_scale*1000));
}

void XFormView::load() {
  QString filename = QFileDialog::getOpenFileName( this, "Choose a file", QDir::currentPath(), "MetaIO/Analyze Files (*.mhd *.mha *.hdr)" );
  if ( !filename.isEmpty() ) {
    pixmap.clear();
    // also clear the imagePoints ...
    imagePoints.clear();

    // First read in the volume ...
    m_volume = new Volume();
    // check format ...
    if (filename.endsWith("mhd"))
      m_volume->InitMhd(filename);
		else if ( filename.endsWith("hdr")) {
			if(! (m_volume->InitAnalyze(filename))) {
				QMessageBox::warning(this, "markTags",
                           "Cannot load Analyze File.\n"
                           "Please check file format.\n\n",
                           "Retry", "Quit", 0, 0, 1);
				return;
			}
		} else {
      QMessageBox::warning(this, "markTags",
                           "Cannot load Volume File.\n"
                           "Please check file format.\n\n",
                           "Retry", "Quit", 0, 0, 1);
      return;
    }

    int bitsAlloc = m_volume->GetVolBitsAllocated();
    int x = m_volume->GetSize(0);
    int y = m_volume->GetSize(1);
    int z = m_volume->GetSize(2);

    m_uiX = x; m_uiY = y;

    // create Pixmaps and change mode ...
    if (bitsAlloc == 8) {
      const unsigned char *pixelData = (unsigned char *)m_volume->GetDataPtr();
      unsigned char *prvu = new unsigned char[4*x*y];
      for (int k=0; k<z; k++) {
        int _max = 0;
        for (int i=0; i<x*y; i++)
          if (_max < pixelData[k*x*y +i])
            _max = pixelData[k*x*y + i];
        float fac = 255.0/_max;

        /* convert the pixel data */

        for (int i=0; i<x*y; i++) {
          prvu[4*i] =   (unsigned char ) (fac*pixelData[k*x*y + i] ) ;
          prvu[4*i+1] = (unsigned char) (fac*pixelData[k*x*y + i]) ;
          prvu[4*i+2] = (unsigned char) (fac*pixelData[k*x*y + i]) ;
          prvu[4*i+3] = 255; // (unsigned char) (fac*pixelData[k*x*y + i]) ;
        }
        // convert to image->pixmap and push back ...
        QImage prv( prvu, x, y, QImage::Format_ARGB32 );
        pixmap.push_back( QPixmap::fromImage(prv) );
      }
      delete [] prvu;
    } else if (bitsAlloc == 16) {
      const unsigned short *pixelData = (unsigned short *)m_volume->GetDataPtr();
      unsigned char *prvu = new unsigned char[4*x*y];
      for (int k=0; k<z; k++) {
        int _max = 0;
        for (int i=0; i<x*y; i++)
          if (_max < pixelData[k*x*y + i])
            _max = pixelData[k*x*y + i];
        float fac = 255.0/_max;

        // std::cout << "Max is " << _max << " and fac is " << fac << std::endl;

        /* convert the pixel data */

        for (int i=0; i<x*y; i++) {
          prvu[4*i] =   (unsigned char ) (fac*pixelData[k*x*y + i] ) ;
          prvu[4*i+1] = (unsigned char) (fac*pixelData[k*x*y + i]) ;
          prvu[4*i+2] = (unsigned char) (fac*pixelData[k*x*y + i]) ;
          prvu[4*i+3] = 255; //(unsigned char) (fac*pixelData[k*x*y + i]) ;
        }
        // convert to image->pixmap and push back ...
        QImage prv( prvu, x, y, QImage::Format_ARGB32 );
        pixmap.push_back( QPixmap::fromImage(prv) );
      }
      delete [] prvu;
    } else {
      QMessageBox::warning(this, "markTags",
                           "Pixel Data type not supported.\n"
                           "Please check file format.\n\n",
                           "Retry", "Quit", 0, 0, 1);
      return;
    }

    emit slicesChanged(z);

    // allocate for imagePoints ...
    for (int i=0; i<z; i++)
      imagePoints.push_back(QPolygonF());

    // adjust scale to display the image properly ...
    QRect r = rect();
    float sx = r.width(); sx /= x;
    float sy = r.height(); sy /= y;
    m_scale = (sx > sy) ? sy : sx;

    m_bHaveVolume = true;
    type = PixmapType; 
    update();
  }
}

void XFormView::reset() {
  QRect r = rect();

  float x = r.width();
  float y = r.height();

  emit rotationChanged(0);
  if (m_bHaveVolume) {
    float sx = r.width(); sx /= m_uiX;
    float sy = r.height(); sy /= m_uiY;
    m_scale = (sx > sy) ? sy : sx;
  } else
    emit scaleChanged(1000);

  ctrlPoints = QPolygonF();
  ctrlPoints << QPointF(x/2, y/2) << QPointF(x/2 + x/5, y/2);
  ctrl->setPoints(ctrlPoints);
  ctrl->firePointChange();


  if ( m_bHaveVolume ) {
    gridPoints = imageToScreen(imagePoints[curr]);
    grid->setPoints(gridPoints);
    // grid->firePointChange();
  }

}

QPolygonF XFormView::imageToScreen(QPolygonF pts) {

  double theta = M_PI*(m_rotation/180.0);
  QPointF center(pixmap[curr].width()/2.0, pixmap[curr].height()/2.0);

  QPointF offset = ctrlPoints.at(0) - center;

  QMatrix trans0(1,0,0,1, offset.x(), offset.y());
  QMatrix trans1(1, 0, 0, 1, center.x(), center.y());
  QMatrix rot(cos(theta), sin(theta), -sin(theta), cos(theta), 0, 0);
  QMatrix scal(m_scale, 0, 0, m_scale, 0, 0);
  QMatrix trans2(1, 0, 0, 1, -center.x(), -center.y());

  QMatrix matrix;
  matrix = trans2*scal*rot*trans1*trans0;

  //std::cout << imagePoints.at(0).x() << " " << imagePoints.at(0).y() << std::endl;

  return matrix.map(pts);
}

QPolygonF XFormView::screenToImage(QPolygonF pts) {
  double theta = M_PI*(m_rotation/180.0);
  QPointF center(pixmap[curr].width()/2.0, pixmap[curr].height()/2.0);

  QPointF offset = ctrlPoints.at(0) - center;

  /*
     QMatrix trans0(1,0,0,1, -offset.x(), -offset.y());
     QMatrix trans1(1, 0, 0, 1, -center.x(), -center.y());
     QMatrix rot(cos(-m_rotation), sin(-m_rotation), -sin(-m_rotation), cos(-m_rotation), 0, 0);
     QMatrix scal(1.0/m_scale, 0, 0, 1.0/m_scale, 0, 0);
     QMatrix trans2(1, 0, 0, 1, center.x(), center.y());


     QMatrix matrix;
     matrix = trans0*trans1*rot*scal*trans2;
     */
  QMatrix trans0(1,0,0,1, offset.x(), offset.y());
  QMatrix trans1(1, 0, 0, 1, center.x(), center.y());
  QMatrix rot(cos(theta), sin(theta), -sin(theta), cos(theta), 0, 0);
  QMatrix scal(m_scale, 0, 0, m_scale, 0, 0);
  QMatrix trans2(1, 0, 0, 1, -center.x(), -center.y());

  QMatrix matrix;
  matrix = trans2*scal*rot*trans1*trans0;

  matrix = matrix.inverted();

  return matrix.map(pts);
}

void XFormView::drawPixmapType(QPainter *painter) {
  QPointF center(pixmap[curr].width()/2.0, pixmap[curr].height()/2.0);
  // Move image to the center ...
  painter->translate(ctrlPoints.at(0) - center);

  // scale and rotate about the center ...
  painter->translate(center);
  painter->rotate(m_rotation);
  painter->scale(m_scale, m_scale);
  painter->translate(-center);

  painter->drawPixmap(QPointF(0, 0), pixmap[curr]);
  painter->setPen(QPen(QColor(255, 0, 0, alpha), 0.25, Qt::SolidLine, Qt::FlatCap, Qt::BevelJoin));
  painter->setBrush(Qt::NoBrush);
  painter->drawRect(QRectF(0, 0, pixmap[curr].width(), pixmap[curr].height()).adjusted(-2, -2, 2, 2));

  // draw grid ...
  if ( m_bHaveVolume && m_bDrawGrid ) {
    painter->setPen(QColor(200,200,0, gridOpacity));
    // Draw grid only if volume is loaded ...
    for (unsigned int i=0; i<m_uiX; i+=gridSpacing) {
      QLineF line(i, 0, i, m_uiY);
      painter->drawLine(line);
    }
    for (unsigned int i=0; i<m_uiY; i+=gridSpacing) {
      QLineF line(0, i, m_uiX, i);
      painter->drawLine(line);
    }
  }
  /*
     std::cout << "--------" << std::endl;
     for (int i=0; i<imagePoints.size(); i++) {
     std::cout << imagePoints.at(i).x() << ", " << imagePoints.at(i).y() << std::endl;
     }
     std::cout << std::endl;
     */
  if ( m_bHaveVolume && m_bDrawPath ) {
    // draw path ...
    // painter->setPen(QColor(60,20,200, gridOpacity));
    painter->setPen(QColor(60,20,200, 255));
    for (int i=0; i< imagePoints[curr].size(); i++) {
      // for each landmark in current frame ...
      // draw FWD ... 
      for (int j=curr+1; j<imagePoints.size(); j++) {
        // stop if this frame doesnt have enough points ...
        if (imagePoints[j].size() <= i) {
          break;
        }
        // else draw the line ....
        QLineF line(imagePoints[j-1].at(i), imagePoints[j].at(i));
        painter->drawLine(line);
      } // end FWD
      // draw BACK ... 
      
      for (int j=(curr-1); j>=0; j--) {
        // stop if this frame doesnt have enough points ...
        if (imagePoints[j].size() <= i) {
          break;
        }
        // else draw the line ....
        QLineF line(imagePoints[j].at(i), imagePoints[j+1].at(i));
        painter->drawLine(line);
      }  // BACK
      
    } // end i
  }
}

void XFormView::drawTextType(QPainter *painter) {
  QFont f("times new roman,utopia");
  f.setStyleStrategy(QFont::ForceOutline);
  f.setPointSize(72);
  f.setStyleHint(QFont::Times);
  painter->setFont(f);

  QFontMetrics fm(f);
  QRectF br(fm.boundingRect(QString("MR Tags")));
  QPointF center(br.center());
  painter->translate(ctrlPoints.at(0) - center);

  painter->translate(center);
  // std::cout << "Rotation is " << m_rotation << std::endl;
  painter->rotate(m_rotation);
  painter->scale(m_scale, m_scale);
  // painter->shear(0, m_shear);
  painter->translate(-center);

  painter->drawText(0, 0, QString("MR Tags"));
  painter->setPen(QPen(QColor(255, 0, 0, alpha), 0.25, Qt::SolidLine, Qt::FlatCap, Qt::BevelJoin));
  painter->setBrush(Qt::NoBrush);
  painter->drawRect(br.adjusted(-1, -1, 1, 1));
}

void XFormView::drawGrid(QPainter *painter) {
  painter->resetMatrix();
  // scale ...
  painter->scale(m_scale, m_scale);

  // set the pen ...
  QPointF offset = QPointF(rect().width()/m_scale, rect().height()/m_scale) - QPointF(m_uiX, m_uiY);

  offset = offset/2;

  painter->setPen(QColor(200,200,0, gridOpacity));

  if ( m_bHaveVolume ) {
    // Draw grid only if volume is loaded ...
    for (int unsigned i=0; i<m_uiX; i+=gridSpacing) {
      QLineF line(i+offset.x(), offset.y(), i+offset.x(), m_uiY+offset.y());
      painter->drawLine(line);
    }
    for (int unsigned i=0; i<m_uiY; i+=gridSpacing) {
      QLineF line(offset.x(), i+offset.y(), m_uiX+offset.x(), i+offset.y());
      painter->drawLine(line);
    }
  }
  // unscale ....
  painter->scale(1.0/m_scale, 1.0/m_scale);

}

void XFormView::drawVectorType(QPainter *painter) {
  QPainterPath path;

  painter->translate(ctrlPoints.at(0) - QPointF(250,250));

  painter->scale(0.77, 0.77);
  painter->translate(98.9154 + 30 , -217.691 - 20);

  QRect br(-55, 275, 500, 590);
  QPoint center = br.center();
  painter->translate(center.x(), center.y());
  painter->rotate(m_rotation);
  painter->scale(m_scale, m_scale);
  // painter->shear(0, m_shear);
  painter->translate(-center.x(), -center.y());
}

void XFormView::copyPrev() {
  // check that copy direction is not reversed.
	if ( m_iLastCopyWasFromNext == 1) {
		int res = QMessageBox::warning ( this, "Changing copy direction", "Last copy was from next. Do you want to continue ?", QMessageBox::Yes, QMessageBox::No);
		if (res == QMessageBox::No)
			return;
	} // if ( m_iLastCopyWasFromNext == 1)

	// can only copy if prev exists
  if (curr) {
		// warn if prev is empty ...
		if ( !imagePoints[curr-1].size() ) {
			int res = QMessageBox::warning ( this, "Empty copy", "Copying from frame with no landmarks. Do you want to continue ?", QMessageBox::Yes, QMessageBox::No);
			if (res == QMessageBox::No)
				return;
		} // if prev is empty 

    imagePoints[curr] = imagePoints[curr-1];
  }
  gridPoints = imageToScreen(imagePoints[curr]);
  grid->setPoints(gridPoints);
	m_iLastCopyWasFromNext = -1;
  update();
}

void XFormView::copyNext() {
	// check that copy direction is not reversed.
	if ( m_iLastCopyWasFromNext == -1) {
		int res = QMessageBox::warning ( this, "Changing copy direction", "Last copy was from previous. Do you want to continue ?", QMessageBox::Yes, QMessageBox::No);
		if (res == QMessageBox::No)
			return;
	}

  // can only copy if next exists
  if (curr < (imagePoints.size()-1) ) {
		// warn if next is empty
		if ( !imagePoints[curr+1].size() ) {
			int res = QMessageBox::warning ( this, "Empty copy", "Copying from frame with no landmarks. Do you want to continue ?", QMessageBox::Yes, QMessageBox::No);
			if (res == QMessageBox::No)
				return;
		} // if next is empty 
    imagePoints[curr] = imagePoints[curr+1];
  }
  gridPoints = imageToScreen(imagePoints[curr]);
  grid->setPoints(gridPoints);
	m_iLastCopyWasFromNext = 1;
  update();
}

void XFormView::loadLM() {
  QString filename = QFileDialog::getOpenFileName( this, "Choose a file", 
													QDir::homePath(), 
													"Landmark Files (*.lm);;Pts(txt) Files (*.pts)" );
  if ( !filename.isEmpty() ) {
    // clear the imagePoints ...
    imagePoints.clear();

	if ( filename.endsWith(".pts") ) {
		// read from txt file
		std::ifstream inFile(filename.toAscii());
		// write out as a txt files ...
		// first the dimensions ... t * n * 2
		int nn, nt;
		float xx, yy;

		inFile >> nt;

		// now the points ...
		for (int i=0; i < nt; i++) {
			inFile >> nn;
			QPolygonF pts;
			for (int j=0; j < nn; j++) {
				inFile >> xx >> yy;
				pts << QPoint(xx, yy);
			}
			imagePoints << pts;
		}
		inFile.close();
	} else {
		// 
		QFile file(filename);
		file.open(QIODevice::ReadOnly);
		QDataStream in(&file);   // we will serialize the data from the file

		in >> imagePoints;

		file.close();
	}
  }  
  update();
}

void XFormView::saveLM() {
  QString filename = QFileDialog::getSaveFileName( this, "Choose a file", 
													QDir::homePath(), 
													"Landmark Files (*.lm);;Pts(txt) Files (*.pts)" );
  if ( !filename.isEmpty() ) {
    
	if ( filename.endsWith(".pts") ) {
		std::ofstream outF(filename.toAscii());
		// write out as a txt files ...
		// first the dimensions ... t * n * 2
		outF << imagePoints.size() << std::endl;

		// now the points ...
		for (int i=0; i < imagePoints.size(); i++) {
			outF << imagePoints.at(i).size() << std::endl;
			for (int j=0; j < imagePoints.at(i).size(); j++) {
				outF << imagePoints.at(i).at(j).x() << " " << imagePoints.at(i).at(j).y() << std::endl; 
			}
		}
		outF.close();
	} else if ( filename.endsWith(".txt") ) {
		// speacial format to be used only with Cardiac for landmarks ...
		unsigned int numpoints = 0;
		for (int i=0; i < imagePoints.size(); i++)
			numpoints += imagePoints.at(i).size();

		std::ofstream outF(filename.toAscii());
		// write out as a txt files ...
		// first the total number of points
		outF << numpoints << " " << imagePoints.size() << std::endl;

		// now write out the points ...
		for (int i=0; i < imagePoints.size(); i++) {
			for (int j=0; j < imagePoints.at(i).size(); j++) {
				outF << imagePoints.at(i).at(j).x() << " " << imagePoints.at(i).at(j).y() << " " << i << std::endl; 
			}
		}
		outF.close();
	
	} else {
		QFile file(filename);
	    file.open(QIODevice::WriteOnly);

		// we will serialize the data from the file
		QDataStream out(&file);   
		out << imagePoints;
		
	    file.close();
	}
  }
}

void XFormWidget::dcm() {
  dcmDialog dlg(this);
  dlg.exec();
}


void XFormWidget::setNumSlices(int n) {
  sliceSlider->setMaximum(n-1);
  sliceSpin->setMaximum(n-1);
}

XFormWidget::XFormWidget(QWidget *parent)
: QWidget(parent) {
  setWindowTitle("Mark Tag Intersections");

  view = new XFormView(this);
  // view->enableOpenGL(true);

  QGroupBox *mainGroup = new QGroupBox(this);
  mainGroup->setFixedWidth(180);
  mainGroup->setTitle("Toolbox");

  QGroupBox *rotateGroup = new QGroupBox(mainGroup);
  rotateGroup->setAttribute(Qt::WA_ContentsPropagated);
  rotateGroup->setTitle("Rotate");
  QSlider *rotateSlider = new QSlider(Qt::Horizontal, rotateGroup);
  rotateSlider->setRange(0, 3600);
  rotateSlider->setValue(0);
  rotateSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

  QGroupBox *scaleGroup = new QGroupBox(mainGroup);
  scaleGroup->setAttribute(Qt::WA_ContentsPropagated);
  scaleGroup->setTitle("Scale");
  QSlider *scaleSlider = new QSlider(Qt::Horizontal, scaleGroup);
  scaleSlider->setRange(1, 12000);
  scaleSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

  QGroupBox *sliceGroup = new QGroupBox(mainGroup);
  sliceGroup->setAttribute(Qt::WA_ContentsPropagated);
  sliceGroup->setTitle(tr("Slice"));
  sliceSlider = new QSlider(Qt::Horizontal, sliceGroup);
  sliceSlider->setRange(0, 30);
  sliceSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
  sliceSpin = new QSpinBox(sliceGroup);
  sliceSpin->setRange(0, 30);
  sliceSpin->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

  QGroupBox *spacingGroup = new QGroupBox(mainGroup);
  spacingGroup->setAttribute(Qt::WA_ContentsPropagated);
  spacingGroup->setTitle(tr("Grid Spacing"));
  QSlider *spacingSlider = new QSlider(Qt::Horizontal, spacingGroup);
  spacingSlider->setRange(1, 20);
  spacingSlider->setValue(4);
  spacingSlider->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
  QSpinBox *spacingSpin = new QSpinBox(spacingGroup);
  spacingSpin->setRange(1, 20);
  spacingSlider->setValue(4);
  spacingSpin->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);

  QGroupBox *copyGroup = new QGroupBox(mainGroup);
  copyGroup->setAttribute(Qt::WA_ContentsPropagated);
  copyGroup->setTitle(tr("Copy Landmarks"));
  QPushButton *copyPrev = new QPushButton(copyGroup);
  copyPrev->setText(tr("Prev"));
  QPushButton *copyNext = new QPushButton(copyGroup);
  copyNext->setText(tr("Next"));

  QGroupBox *typeGroup = new QGroupBox(mainGroup);
  typeGroup->setAttribute(Qt::WA_ContentsPropagated);
  typeGroup->setTitle(tr("View"));
  QCheckBox *hoverView = new QCheckBox(typeGroup);
  hoverView->setCheckState(Qt::Checked);
  QCheckBox *tagView = new QCheckBox(typeGroup);
  tagView->setCheckState(Qt::Unchecked);
  hoverView->setText(tr("Control"));
  tagView->setText(tr("Grid"));
  QCheckBox *pathView = new QCheckBox(typeGroup);
  pathView->setCheckState(Qt::Unchecked);
  pathView->setText(tr("Path"));

  QGroupBox *dataGroup = new QGroupBox(mainGroup);
  dataGroup->setAttribute(Qt::WA_ContentsPropagated);
  dataGroup->setTitle("Data");

  QPushButton *dcmButton = new QPushButton(dataGroup);
  dcmButton->setText("Import DICOM");
  QPushButton *loadButton = new QPushButton(dataGroup);
  loadButton->setText("Load Volume");
  QPushButton *saveLMButton = new QPushButton(dataGroup);
  saveLMButton->setText("Save LM");
  QPushButton *loadLMButton = new QPushButton(dataGroup);
  loadLMButton->setText("Load LM");

  QPushButton *resetButton = new QPushButton(mainGroup);
  resetButton->setText("Reset Transform");

  // LAYOUTS ...

  QHBoxLayout *viewLayout = new QHBoxLayout(this);
  viewLayout->addWidget(view);
  viewLayout->addWidget(mainGroup);

  QHBoxLayout *rotateGroupLayout = new QHBoxLayout(rotateGroup);
  rotateGroupLayout->addWidget(rotateSlider);

  QHBoxLayout *scaleGroupLayout = new QHBoxLayout(scaleGroup);
  scaleGroupLayout->addWidget(scaleSlider);

  QHBoxLayout *sliceGroupLayout = new QHBoxLayout(sliceGroup);
  sliceGroupLayout->addWidget(sliceSlider);
  sliceGroupLayout->addWidget(sliceSpin);

  QHBoxLayout *spacingGroupLayout = new QHBoxLayout(spacingGroup);
  spacingGroupLayout->addWidget(spacingSlider);
  spacingGroupLayout->addWidget(spacingSpin);

  QHBoxLayout *copyLayout = new QHBoxLayout(copyGroup);
  copyLayout->addWidget(copyPrev);
  copyLayout->addWidget(copyNext);

  QGridLayout *typeGroupLayout = new QGridLayout(typeGroup);
  typeGroupLayout->addWidget(hoverView,1,1);
  typeGroupLayout->addWidget(tagView,2,1);
  typeGroupLayout->addWidget(pathView,1,2);

  QVBoxLayout *dataGroupLayout = new QVBoxLayout(dataGroup);
  dataGroupLayout->addWidget(dcmButton);
  dataGroupLayout->addWidget(loadButton);
  dataGroupLayout->addWidget(loadLMButton);
  dataGroupLayout->addWidget(saveLMButton);

  QVBoxLayout *mainGroupLayout = new QVBoxLayout(mainGroup);
  mainGroupLayout->addWidget(rotateGroup);
  mainGroupLayout->addWidget(scaleGroup);
  mainGroupLayout->addWidget(resetButton);

  mainGroupLayout->addWidget(sliceGroup);
  mainGroupLayout->addWidget(spacingGroup);
  mainGroupLayout->addStretch(1);
  mainGroupLayout->addWidget(copyGroup);
  mainGroupLayout->addStretch(1);
  mainGroupLayout->addWidget(typeGroup);
  mainGroupLayout->addWidget(dataGroup);

  // SIGNALS-SLOTS

  connect(rotateSlider, SIGNAL(valueChanged(int)), view, SLOT(changeRotation(int)));
  connect(scaleSlider, SIGNAL(valueChanged(int)), view, SLOT(changeScale(int)));
  connect(sliceSlider, SIGNAL(valueChanged(int)), view, SLOT(changeSlice(int)));
  connect(sliceSlider, SIGNAL(valueChanged(int)), sliceSpin, SLOT(setValue(int)));
  connect(sliceSpin, SIGNAL(valueChanged(int)), sliceSlider, SLOT(setValue(int)));
  connect(spacingSlider, SIGNAL(valueChanged(int)), view, SLOT(changeSpacing(int)));
  connect(spacingSlider, SIGNAL(valueChanged(int)), spacingSpin, SLOT(setValue(int)));
  connect(spacingSpin, SIGNAL(valueChanged(int)), spacingSlider, SLOT(setValue(int)));

  connect(view, SIGNAL(rotationChanged(int)), rotateSlider, SLOT(setValue(int)));
  connect(view, SIGNAL(scaleChanged(int)), scaleSlider, SLOT(setValue(int)));

  connect(copyPrev, SIGNAL(clicked()), view, SLOT(copyPrev()));
  connect(copyNext, SIGNAL(clicked()), view, SLOT(copyNext()));

  connect(resetButton, SIGNAL(clicked()), view, SLOT(reset()));
  connect(view, SIGNAL(descriptionEnabledChanged(bool)), view->hoverPoints(), SLOT(setDisabled(bool)));
  
  connect(dcmButton, SIGNAL(clicked()), this, SLOT(dcm()));
  connect(loadButton, SIGNAL(clicked()), view, SLOT(load()));
  connect(loadLMButton, SIGNAL(clicked()), view, SLOT(loadLM()));
  connect(saveLMButton, SIGNAL(clicked()), view, SLOT(saveLM()));

  connect(view, SIGNAL(slicesChanged(int)), this, SLOT(setNumSlices(int)));
  connect(hoverView, SIGNAL(stateChanged(int)), view, SLOT(enableHover(int)));
  connect(tagView, SIGNAL(stateChanged(int)), view, SLOT(enableGrid(int)));
  connect(pathView, SIGNAL(stateChanged(int)), view, SLOT(enablePath(int)));

  // defaults
  emit resetButton->animateClick();
}
