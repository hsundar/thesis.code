#include "viewer3d.h"
#include "Volume.h"

#include <fstream>
#include <iostream>

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoBlinker.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTexture3.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/manips/SoClipPlaneManip.h>

#include <VolumeViz/nodes/SoTransferFunction.h>
#include <VolumeViz/nodes/SoVolumeData.h>
#include <VolumeViz/nodes/SoVolumeRender.h>
#include <VolumeViz/nodes/SoVolumeRendering.h>

// #include "colorswatch.h"
#include "tooldock.h"
#include "mainwindow.h"
#include "MyExaminerViewer.h"
#include "BoundingBox.h"
#include "MarchingCubes.h"
#include "FiberTracker.h"
#include "Volume.h"
#include "Point.h"

// vec field display ..
#include "SFVectorField.h"
#include "SoVectorFieldViz.h"

viewer3d::viewer3d(QWidget* parent, Qt::WFlags f):
QWidget(parent, f), m_iDepth(8), m_bHaveVolume(false), m_bHaveOctree(false), m_bHaveField(false) {

	// QFrame *frame = new QFrame(this);
	// Initialize the libraries ...
	SoQt::init(this);
	SoVolumeRendering::init();
	// this is deprecated. move to rendering a simple cube instead.
	BoundingBox::initClass();
	// Initialize new classes
	SFScalarField::initClass();
	SFVectorField::initClass();
	MarchingCubes::initClass();
	SoVectorFieldViz::initClass();

	m_soBoundingBox = NULL;
	// m_soVolumeData = NULL;
	m_soTF = NULL;

	// Create the empty root scenegraph
	m_sceneGraph = new SoSeparator;
	m_sceneGraph->ref();

	m_soVolSwitch = new SoBlinker;
	m_soVolSwitch->on = FALSE; 

	m_volGroup = new SoSeparator;
	m_volGroup->addChild(m_soVolSwitch);

	m_geoGroup = new SoSeparator;
	m_geoGroup->ref();
	SoBaseColor * col = new SoBaseColor;
	col->rgb = SbColor(1, 1, 0);
	m_geoGroup->addChild(col);

	m_sceneGraph->addChild(m_volGroup);
	m_sceneGraph->addChild(m_geoGroup);
	m_soOctSwitch = new SoSwitch;
	m_soOctSwitch->ref();
	m_geoGroup->addChild(m_soOctSwitch);

	// m_volume = NULL;

	// Initialize the colormap
	m_ucColors = new unsigned char[4*256]; // @bug hardcoded for now. 8-bit RGBA colormap.
	for (unsigned int i=0; i<50; i++)
		m_ucColors[i*4+3]=0;
	for (unsigned int i=50; i<256; i++)
		m_ucColors[i*4+3]=50;

	m_tracker = new FiberTracker();

	// Now lets add the examiner viewer ... 
	m_eView = new MyExaminerViewer (this); //MyExaminerViewer(this);

	// m_eViewer->setTitle("Hari's Volume Viewer");
	m_eView->setBackgroundColor(SbColor(0.176f, 0.234f, 0.477f));
	m_eView->setSceneGraph(m_sceneGraph);
	m_eView->setFeedbackVisibility(true);
	// eView->setDrawStyle(SoQtViewer::STILL, SoQtViewer::VIEW_WIREFRAME_OVERLAY);
	m_eView->setDrawStyle(SoQtViewer::INTERACTIVE, SoQtViewer::VIEW_LOW_COMPLEXITY);
	m_eView->setTransparencyType(SoGLRenderAction::BLEND);
	m_eView->setStereoOffset(200.0f);

	m_defViz=NULL;

	m_eView->show();
}

void viewer3d::loadVolume(const QString &filename) {
	// First read in the volume ...
	Volume *_volume = new Volume();
	// check format ...
	if (filename.endsWith("mhd"))
		_volume->InitMhd(filename);
	else if ( filename.endsWith("hdr"))
		_volume->InitAnalyze(filename);
	else {
		QMessageBox::warning(this, "Cardiac Viewer",
			"Cannot load Volume File.\n"
			"Please check file format.\n\n",
			"Retry", "Quit", 0, 0, 1);
		return;
	}

	//std::cout << "Finished loading Volume" << std::endl;
	//~~~~~~~~~~~~~~~~~~~
	// compute histogram ... and set it ...
	// int histSize = 1 << m_volume->GetVolBitsStored();
	// unsigned int *hist = new unsigned int[histSize];

	// m_volume->ComputeHistogram(hist, histSize);

	// fix this to set the new TFE 
	// m_tfe->setHistogram(hist, histSize);
	// m_tfe->update();

	// delete [] hist;

	//~~~~~~~~~~~~~~~~~~~~
	SbVec3s dim = SbVec3s(_volume->GetSize(0), _volume->GetSize(1), _volume->GetSize(2));
	//_volSize = Point(_volume->GetSize(0), _volume->GetSize(1), _volume->GetSize(2));
	//_volSpacing = Point(_volume->GetVoxelSize(0), _volume->GetVoxelSize(1), _volume->GetVoxelSize(2));

	SoSFVec3f _mn;
	_mn.setValue( - _volume->GetVoxelSize(0)*_volume->GetSize(0)/2, - _volume->GetVoxelSize(1)* _volume->GetSize(1)/2, - _volume->GetVoxelSize(2)* _volume->GetSize(2)/2 );
	SoSFVec3f _mx;
	_mx.setValue( _volume->GetVoxelSize(0)* _volume->GetSize(0)/2, _volume->GetVoxelSize(1)* _volume->GetSize(1)/2, _volume->GetVoxelSize(2)*_volume->GetSize(2)/2 );

	SbBox3f bbox(- _volume->GetVoxelSize(0)* _volume->GetSize(0)/2, - _volume->GetVoxelSize(1)* _volume->GetSize(1)/2, \
		- _volume->GetVoxelSize(2)* _volume->GetSize(2)/2, _volume->GetVoxelSize(0)* _volume->GetSize(0)/2, \
		_volume->GetVoxelSize(1)* _volume->GetSize(1)/2,  _volume->GetVoxelSize(2)* _volume->GetSize(2)/2 ); 

	// set the labels ... 
	// Maybe set this within the Inventor window ?
	// Set this later ...
	// @bug -- not the right way to set ... send signal instead.
	toolDock *dock = ((MainWindow *)parentWidget())->getToolDock();
	QString dimStr = QString::number( _volume->GetSize(0)) + QString("x") + QString::number( _volume->GetSize(1)) + QString("x") + QString::number( _volume->GetSize(2)) + QString(" <b>voxels</b>");
	// std::cout << qPrintable(dimStr) << std::endl;
	dock->setVolDims( dimStr );
	dock->setVolSpacing( QString::number( _volume->GetVoxelSize(0), 'f', 2) + QString("x") + QString::number( _volume->GetVoxelSize(1), 'f', 2) + QString("x") + QString::number( _volume->GetVoxelSize(2), 'f', 2) + QString(" <b>mm</b>"));

	// Now create a SoVolumeNode ..
	// Add SoVolumeData to scene graph
	SoVolumeData *_soVolumeData = new SoVolumeData();
	_soVolumeData->setVolumeData(dim, _volume->GetDataPtr(), SoVolumeData::UNSIGNED_BYTE); 
	_soVolumeData->setVolumeSize(bbox);

	// add to the lists ...
	m_volume.append(_volume);
	m_soVolumeData.append(_soVolumeData);

	// the scenegraph ...
	SoSeparator *volFrameSep = new SoSeparator;

	// add the isosurface ,,,// test isoSurface ...
	float sx1, sy1, sz1;
	float tx1, ty1, tz1;

	bbox.getSize( sx1, sy1, sz1 );
	bbox.getOrigin( tx1, ty1, tz1 );

	
	SoSeparator *isoSep = new SoSeparator;
	SoTransform *isoTrans = new SoTransform;
	isoTrans->scaleFactor.setValue(sx1, sy1, sz1);
	isoTrans->translation.setValue(tx1, ty1, tz1);
	isoSep->addChild(isoTrans);
	isoSep->addChild(addIsosurface(_volume));
	volFrameSep->addChild(isoSep);
	
	// m_soVolSwitch->addChild(_soVolumeData);

	// The actual bounding box ...
	if (m_soBoundingBox == NULL) {
		m_soBoundingBox = new BoundingBox();
		m_soBoundingBox->textOn = false;
		m_soBoundingBox->setName("BoundingBox");

		m_soBBoxSwitch = new SoSwitch;
		m_soBBoxSwitch->whichChild = SO_SWITCH_ALL;
		m_soBBoxSwitch->setName("BoundingBoxSwitch");

		m_soBBoxSwitch->addChild(m_soBoundingBox);
		m_geoGroup->addChild(m_soBBoxSwitch);
	}

	m_soBoundingBox->minBounds = _mn;
	m_soBoundingBox->maxBounds = _mx;

	// change Octree transform if already loaded ...
	if (m_bHaveOctree) {
		SoTransform *o2v = (SoTransform *)SoNode::getByName("oct2Vol");
		float sx, sy, sz;
		float tx, ty, tz;

		float lng = (float)(1 << m_iDepth);
		bbox.getSize( sx, sy, sz );
		sx /= lng; sy /= lng; sz /= lng;
		bbox.getOrigin( tx, ty, tz );

		o2v->scaleFactor.setValue(sx, sy, sz);
		o2v->translation.setValue(tx, ty, tz);

		SoDrawStyle *style = (SoDrawStyle *)SoNode::getByName("octStyle");
		style->style.setValue(SoDrawStyle::LINES);
	}

	if ( m_soVolumeData.size() < 2 ) {
		// Add TransferFunction (color map) to scene graph
		m_soTF = new SoTransferFunction();
		m_soTF->setName("TransferFuction");

		for (int i=0;i<256*4;++i) {
			m_soTF->colorMap.set1Value(i, m_ucColors[i]/256.0f);
		}

		// Add VolumeRender to scene graph
		m_soVolRen = new SoVolumeRender();

		// m_volGroup->addChild(m_soVolRen); 
		// m_volGroup->addChild(m_soTF);  
	} else {
		// m_soVolSwitch->on = TRUE;
		m_soVolSwitch->speed = 0.2;
	}

	volFrameSep->addChild(_soVolumeData);
	volFrameSep->addChild(m_soTF);
	volFrameSep->addChild(m_soVolRen);

	m_soVolSwitch->addChild(volFrameSep);

	//SoGetBoundingBoxAction ba(m_eView->getViewportRegion());
	//ba.apply(m_sceneGraph);

	//SbBox3f box = ba.getBoundingBox();

	// @bug should not be added all the time ...

	m_soClipSwitch = new SoSwitch;
	m_soClipSwitch->whichChild = SO_SWITCH_NONE;
	m_soClipSwitch->setName("ClipPlaneSwitch");
	m_soClipPlane = new SoClipPlaneManip;
	m_soClipPlane->setValue(bbox, SbVec3f(0.0f, 0.0f, -1.0f), 1.1f);
	m_soClipSwitch->addChild(m_soClipPlane); 
	m_volGroup->insertChild(m_soClipSwitch, 0);
	// m_geoGroup->insertChild(m_soClipSwitch, 0);

	m_bHaveVolume = true;
	m_soVolSwitch->whichChild = 0;

	dock->showTab(0);
	m_eView->viewAll();
	m_eView->saveHomePosition();
	// std::cout << "Leaving loadVolume" << std::endl;
}

void viewer3d::colormapChanged() {
	if (m_bHaveVolume) {
		m_soTF->predefColorMap.setValue(SoTransferFunction::NONE);
		m_soTF->colorMapType.setValue(SoTransferFunction::RGBA);

		// qDebug ( "COLORMAP %d, %d, %d and alpha is %d", m_ucColors[0], m_ucColors[1], m_ucColors[2], m_ucColors[3]);
		for (int i=0;i<256*4;++i)
			m_soTF->colorMap.set1Value(i, m_ucColors[i]/256.0f);

		// m_soTF->reMap(1, 256);
		m_soTF->reMap(1, 65535);
	}
}

void viewer3d::loadOctree( const QString &filename) {
	int dim, depth, n;

	std::ifstream in (filename.toAscii());
	in >> dim;
	in >> depth;
	in >> n;

	m_iDepth = depth;

	// ((MainWindow*)parent())->getDock()->setMaxDepth(depth);
	toolDock *dock = ((MainWindow *)parentWidget())->getToolDock();
	dock->setMaxDepth(depth);
	// spinHigh->setMaxValue(depth);
	// spinLow->setMaxValue(depth);

	float lng = (float)(1<<(depth));

	m_eView->setMaxDepth(depth);

	// center the octree ...
	SoTransform *oct2Vol = new SoTransform;
	oct2Vol->setName("oct2Vol");

	SoDrawStyle * style = new SoDrawStyle;
	style->setName("octStyle");

	SoSFVec3f _mn;
	SoSFVec3f _mx;

	if (m_bHaveVolume) {
		float sx, sy, sz;
		float tx, ty, tz;

		SbBox3f bbox = m_soVolumeData[0]->getVolumeSize();

		bbox.getSize( sx, sy, sz );
		sx /= lng; sy /= lng; sz /= lng;
		bbox.getOrigin( tx, ty, tz );

		_mn = bbox.getMin();
		_mx = bbox.getMax();

		oct2Vol->scaleFactor.setValue(sx, sy, sz);
		oct2Vol->translation.setValue(tx, ty, tz);

		style->style.setValue(SoDrawStyle::LINES);
	} else {
		// oct2Vol->translation.setValue(-lng/2, -lng/2, -lng/2);
		_mn.setValue(-lng/2, -lng/2, -lng/2);
		_mx.setValue(lng/2, lng/2, lng/2);
		_mn.setValue(0, 0, 0);
		_mx.setValue(lng, lng, lng);
		// style->style.setValue(SoDrawStyle::FILLED);
		style->style.setValue(SoDrawStyle::LINES);
	}
	srand(1713);
	for (int i=0; i<depth+1; i++) {
		// std::cout << "i = " << i << std::endl;
		SoSeparator *sep = new SoSeparator;
		SoSwitch *sw = new SoSwitch();
		sep->addChild(style);
		sep->addChild(oct2Vol);

		// add name ...
		char str[32];
		sprintf(str, "depth%d", i);
		sep->setName(str);
		sprintf(str, "switch%d", i);
		sw->setName(str);
		// level[i]->ref();
		m_soOctSwitch->addChild(sw);
		// m_geoGroup->addChild(sw);
		sw->addChild(sep);
		sw->whichChild = SO_SWITCH_ALL;
		SoMaterial *mat = new SoMaterial;
		float val = (float)(i)/(depth+1);
		/*
		if (i < 3) val = 0;
		if (i==4) val = 0.2;
		if (i==5) val = 0.25;
		if (i==6) val = 0.4;
		*/
		// std::cout << "Val is " << val << std::endl;
		mat->transparency.setValue(1.0-val);
		float r = ((float)rand()/RAND_MAX);
		float g = ((float)rand()/RAND_MAX);
		float b = ((float)rand()/RAND_MAX);
		mat->ambientColor.setValue(r, g, b);
		mat->diffuseColor.setValue(r, g, b);
		mat->specularColor.setValue(r, g, b); // 0.9, 0.9, 0.5);
		mat->emissiveColor.setValue(r, g, b);
		mat->shininess = 0.3;
		sep->addChild(mat);
	} 

	SoBaseColor * col = new SoBaseColor;
	col->rgb = SbColor(1, 0, 0);
	m_sceneGraph->addChild(col);

	// The actual bounding box ...
	if (m_soBoundingBox == NULL) {
		m_soBoundingBox = new BoundingBox();
		m_soBoundingBox->textOn = false;
		m_soBoundingBox->setName("BoundingBox");

		m_soBBoxSwitch = new SoSwitch;
		m_soBBoxSwitch->whichChild = SO_SWITCH_ALL;
		m_soBBoxSwitch->setName("BoundingBoxSwitch");

		m_soBBoxSwitch->addChild(m_soBoundingBox);
		m_geoGroup->addChild(m_soBBoxSwitch);
	}

	m_soBoundingBox->minBounds = _mn;
	m_soBoundingBox->maxBounds = _mx;


	int x,y,z,d;

	for (int i=0; i<n; i++) {
		// read in a point ...
		in >> x >> y >> z >> d;
		// std::cout << i << ":  " << x << " " << y << " "  << z  << " " << d << std::endl;

		if ( d > 6)
			continue;

		//  if ( ( z < 80 ) ) // || (z > 90) )
		//    continue;

		char str[32];
		sprintf(str, "depth%d", d);
		SoSeparator * level = (SoSeparator *)SoNode::getByName(str);

		// add the appropriate cube .. 
		int sz = 1 << (depth - d);
		// std::cout << "size is " << sz << std::endl;
		SoSeparator *sep = new SoSeparator;
		// std::cout << "adding separator under level ... " << std::endl;
		level->addChild(sep);
		float s1,s2,s3; s1=sz;
		// std::cout << "setting scale" << std::endl;
		SoTransform *trans = new SoTransform;
		//trans->scaleFactor.setValue(s1, s1, s1);
		s1=(float)sz/2 +x; s2=(float)sz/2 + y; s3=(float)sz/2 + z;
		// std::cout << "setting translation" << std::endl;
		trans->translation.setValue(s1, s2, s3);
		sep->addChild(trans);
		s1=sz;
		// std::cout << "adding cube under separator under level ... " << std::endl;
		SoCube *cube = new SoCube;
		cube->width = sz; cube->height = sz; cube->depth = sz;
		// std::cout << "Added cube: " << sz << " " << s1 << " " << s2 << " " << s3 << std::endl;
		sep->addChild(cube);
	}

	// m_geoGroup->addChild(bb);

	m_bHaveOctree = true;
	m_soOctSwitch->whichChild = SO_SWITCH_ALL;

	dock->showTab(1);

	m_eView->viewAll();
	m_eView->saveHomePosition();
}

void viewer3d::toggleClipPlane() {
	if ( !m_bHaveVolume )
		return;
	int val = m_soClipSwitch->whichChild.getValue();
	m_soClipSwitch->whichChild = (val == SO_SWITCH_NONE) ? SO_SWITCH_ALL : SO_SWITCH_NONE;
}

void viewer3d::toggleBoundingBox() {
	if ( !m_bHaveVolume && !m_bHaveOctree )
		return;
	int val = m_soBBoxSwitch->whichChild.getValue();
	m_soBBoxSwitch->whichChild = (val == SO_SWITCH_NONE) ? SO_SWITCH_ALL : SO_SWITCH_NONE;
}

void viewer3d::toggleVolume() {
	int val = m_soVolSwitch->whichChild.getValue();
	m_soVolSwitch->whichChild = (val == SO_SWITCH_NONE) ? SO_SWITCH_ALL : SO_SWITCH_NONE;
}

void viewer3d::playSeq() {
	controlSequencer(PLAY);
}

void viewer3d::pauseSeq() {
	controlSequencer(PAUSE);
}

void viewer3d::rewSeq() {
	controlSequencer(REW);
}

void viewer3d::fwdSeq() {
	controlSequencer(FWD);
}

void viewer3d::controlSequencer(playType mode) {
	int val;
	switch (mode) {
	case PLAY:
		m_soVolSwitch->on = TRUE;
		break;
	case PAUSE:
		m_soVolSwitch->on = FALSE;
		break;
	case REW:
		val = m_soVolSwitch->whichChild.getValue() -1;
		if (val < 0) val = m_soVolumeData.size() -1;
		m_soVolSwitch->whichChild = val;
		break;
	case FWD:
		val = m_soVolSwitch->whichChild.getValue() + 1;
		if (val >= m_soVolumeData.size()) val = 0;
		m_soVolSwitch->whichChild = val;
		break;
	}
}

void viewer3d::toggleOctree() {
	int val = m_soOctSwitch->whichChild.getValue();
	m_soOctSwitch->whichChild = (val == SO_SWITCH_NONE) ? SO_SWITCH_ALL : SO_SWITCH_NONE;
}

void viewer3d::setOctThreshold(int l, int h) {
	if (! m_bHaveOctree )
		return;

	for (int i=0; i<m_iDepth+1; i++) {
		// get the right switch ...
		char str[32];
		sprintf(str, "switch%d", i);
		SoSwitch * sw = (SoSwitch *)SoNode::getByName(str);

		if ( (i >= l) && (i <= h) )
			sw->whichChild = SO_SWITCH_ALL;
		else
			sw->whichChild = SO_SWITCH_NONE;
	}
}

void viewer3d::generateOctree() {
	if (m_bHaveVolume)
		m_volume[0]->ComputeLinearOctree();

	// now render it ...
}

void viewer3d::redrawDeformationField() {
	if (m_defField.size()) {
		// get params from the toolDock
		toolDock *dock = ((MainWindow *)parentWidget())->getToolDock();

		int idx = dock->defFieldCombo->currentIndex();
		int dim[3]; dim[0] = m_defField[idx]->getSize(0); dim[1] = m_defField[idx]->getSize(1); dim[2] = m_defField[idx]->getSize(2);

		m_defViz->data.setValue(dim, m_defField[idx]->getDataPtr(), false);

		int plane = dock->defVisPlaneCombo->currentIndex();
		int density = dock->defVisSparsitySlider->value();
		dock->defVisSliceSlider->setRange(0,dim[plane]-1);
		dock->defVisSliceSpin->setRange(0,dim[plane]-1);

		m_defViz->sparsity = 1./density;
		m_defViz->drawMode = dock->defVisModeCombo->currentIndex();
		m_defViz->colorMode = dock->defVisColorCombo->currentIndex();
		m_defViz->gridPlane = dock->defVisPlaneCombo->currentIndex();
		m_defViz->sliceNum = dock->defVisSliceSlider->value();
		return;
	}

}

void viewer3d::loadField (const QString &fname) {
	Field<float, 3>* _fld = new Field<float, 3>();
	m_defField.append(_fld);
	int idx = m_defField.size() -1;
	// bool flag = m_pPD->init(fname.toStdString());
	m_defField[idx]->init(fname.toStdString());  
	
	// get params from the toolDock
	toolDock *dock = ((MainWindow *)parentWidget())->getToolDock();
	
	
	int dim[3]; dim[0] = m_defField[idx]->getSize(0); dim[1] = m_defField[idx]->getSize(1); dim[2] = m_defField[idx]->getSize(2);
	int dim2[3]; dim2[0] = m_volume[0]->GetSize(0); dim2[1] = m_volume[0]->GetSize(1); dim2[2] = m_volume[0]->GetSize(2);

	
	int plane = dock->defVisPlaneCombo->currentIndex();
	int density = dock->defVisSparsitySlider->value();
	dock->defVisSliceSlider->setRange(0,dim[plane]-1);
	dock->defVisSliceSpin->setRange(0,dim[plane]-1);

	if (m_bHaveField) {
		dock->defFieldCombo->addItem(QString::number(idx));
		return;
	}

	// create the visualizer ...
	m_defViz = new SoVectorFieldViz;

	SoSeparator *root = new SoSeparator;

	// center the octree ...
	SoTransform *pts2Vol = new SoTransform;
	pts2Vol->setName("pts2Vol");

	SoSFVec3f _mn;
	SoSFVec3f _mx;

	float* buff = NULL;

	if (m_bHaveVolume) {
		float sx, sy, sz;
		float tx, ty, tz;

		SbBox3f bbox = m_soVolumeData[0]->getVolumeSize();

		bbox.getSize( sx, sy, sz );
		sx /= dim[0]; sy /= dim[1]; sz /= dim[2];
		bbox.getOrigin( tx, ty, tz );

		_mn = bbox.getMin();
		_mx = bbox.getMax();

		pts2Vol->scaleFactor.setValue(sx, sy, sz);
		pts2Vol->translation.setValue(tx, ty, tz);

		// create the image buffer ...
		buff = new float[dim[0]*dim[1]*dim[2]];
		unsigned char* vol  = (unsigned char*)m_volume[0]->GetDataPtr();
		
		for (unsigned int i=0; i<dim[0]*dim[1]*dim[2]; i++)
			buff[i] = 0.0; 
		for (unsigned int k=0; k<dim2[2]; k++)
			for (unsigned int j=0; j<dim2[1]; j++)
				for (unsigned int i=0; i<dim2[0]; i++) {
					buff[(k*dim[1]+j)*dim[0]+i] = static_cast<float>(vol[(k*dim2[1]+j)*dim2[0]+i]);
				}
	} 

	root->addChild(pts2Vol);

	/*  SoTransform *defTrans = new SoTransform;
	defTrans->scaleFactor.setValue(sx1, sy1, sz1);
	defTrans->translation.setValue(tx1, ty1, tz1);
	defSep->addChild(isoTrans);
	*/

	// Set the field data
	m_defViz->data.setValue(dim, m_defField[idx]->getDataPtr(), false);
	// set the image data -- copy the image info.
	m_defViz->image.setValue(dim, buff, true);
	// set the fiber orientations 

	m_defViz->sparsity = 1./density;
	m_defViz->drawMode = dock->defVisModeCombo->currentIndex();
	m_defViz->gridPlane = dock->defVisPlaneCombo->currentIndex();
	m_defViz->sliceNum = dock->defVisSliceSlider->value();

	// Create a VertexProperty because they are more efficient the
	// SoCoordinate3 nodes.
	SoVertexProperty *vprop = new SoVertexProperty();
	vprop->ref();
	vprop->materialBinding = SoVertexProperty::PER_VERTEX;
	vprop->vertex.connectFrom(&m_defViz->points);
	vprop->orderedRGBA.connectFrom(&m_defViz->colors);

	// connect the indexes from MarchingCubes
	SoLineSet *lineSet = new SoLineSet();
	lineSet->vertexProperty = vprop;
	lineSet->numVertices.connectFrom(&m_defViz->lineSet);
	root->addChild(lineSet);

	m_sceneGraph->addChild(root);
	// m_soVolSwitch->addChild(root);

	dock->showTab(1);
	
	m_bHaveField = true;
	dock->defFieldCombo->addItem(QString::number(idx));

	m_eView->viewAll();
	m_eView->saveHomePosition();

	return;
}

SoSeparator* viewer3d::addIsosurface(Volume* vol) {
	SoSeparator *root = new SoSeparator();

	// iso material ...
	SoMaterial *mat = new SoMaterial;

	mat->transparency.setValue(0.0);
	float r = 0.8;
	float g = 0.2;
	float b = 0.2;
	mat->ambientColor.setValue(r, g, b);
	mat->diffuseColor.setValue(r, g, b);
	mat->specularColor.setValue(r, g, b); 
	mat->emissiveColor.setValue(r, g, b);
	mat->shininess = 0.3;

	//  SoBaseColor * col = new SoBaseColor;
	//  col->rgb = SbColor(1, 0, 0);
	root->addChild(mat);

	// create MarchingCubes engine and connect fields
	MarchingCubes *mcubes = new MarchingCubes();
	mcubes->ref();

	// set the data for mcubes ...
	// for now convert it to float ...
	// @bug  ... chaneg Marching cubes to allow this later ...
	int dim[3]; dim[0] = vol->GetSize(0); dim[1] = vol->GetSize(1); dim[2] = vol->GetSize(2);
	long sz = dim[0]* dim[1]* dim[2];
	float *data = new float[sz];
	unsigned char *volData = (unsigned char *) vol->GetDataPtr();

	for (int i=0; i<sz; i++)
		data[i] = (float)(volData[i]);

	// qDebug("Setting data for marching cubes ...");
	mcubes->data.setValue(dim, data, false);
	// qDebug("Finished setting data for marching cubes");

	mcubes->isoValue = 80.;
	// now to generate ...
	// MarchingCubes gives the triangles with the vertices ordered clockwise
	SoShapeHints *hints = new SoShapeHints();
	hints->vertexOrdering.setValue(SoShapeHints::CLOCKWISE);
	root->addChild(hints);

	// Create a VertexProperty because they are more efficient the
	// SoCoordinate3 nodes.
	SoVertexProperty *vprop = new SoVertexProperty();
	vprop->ref();
	vprop->vertex.connectFrom(&mcubes->points);
	vprop->normal.connectFrom(&mcubes->normals);

	// connect the indexes from MarchingCubes
	SoIndexedFaceSet *faceSet = new SoIndexedFaceSet();
	faceSet->vertexProperty = vprop;
	faceSet->coordIndex.connectFrom(&mcubes->indexes);
	root->addChild(faceSet);

	// add to the list ...
	m_isoMaterial.append(mat);
	m_marchers.append(mcubes);

	return root;
}

void viewer3d::changeIsoLevel(int lev) {
	if (!m_bHaveVolume)
		return;
	// change for current volume ... and update later ...
	int val = m_soVolSwitch->whichChild.getValue();

	m_marchers[val]->isoValue = (float)lev;
	/*
	QListIterator<MarchingCubes*> i(m_marchers);
	while (i.hasNext())
	i.next()->isoValue = (float)lev;
	*/
}

void viewer3d::updateIsoAll(int lev) {
	if (!m_bHaveVolume)
		return;
	QListIterator<MarchingCubes*> i(m_marchers);

	int cnt=0;
	QProgressDialog progress("Computing  Isosurfaces ...", "Abort", 0, m_marchers.size(), this);
	progress.setValue(cnt);
	while (i.hasNext()) {
		progress.setValue(cnt++);
		qApp->processEvents();
		i.next()->isoValue = (float)lev;

		if (progress.wasCanceled())
			break;
	}
	progress.setValue(m_marchers.size());
}

// for now does a complete clean up ...
void viewer3d::unloadVolume() {


}

void viewer3d::loadFA(const QString &fname) {
	m_pFA = new Field<float, 1>();
	// bool flag = m_pFA->init(fname.toStdString());
	m_pFA->init(fname.toStdString());
	/*	if (flag)
	statusBar()->showMessage(tr("Successfully loaded Fractional Anisotropy Image"), 2000);
	else
	statusBar()->showMessage(tr("Could not load Fractional Anisotropy Image"), 2000);
	*/
}

void viewer3d::loadPD(const QString &fname) {
	m_pPD = new Field<float, 3>();
	// bool flag = m_pPD->init(fname.toStdString());
	m_pPD->init(fname.toStdString());
	/*	if (flag)
	statusBar()->showMessage(tr("Successfully loaded Principal Directions Image"), 2000);
	else
	statusBar()->showMessage(tr("Could not load Principal Directions Image"), 2000);
	*/
}

void viewer3d::fiberTrack(int density, int color) {
	FiberTracker::Color clr = FiberTracker::Normal;
	switch (color) {
	case 0: clr = FiberTracker::Gray; break;
	case 1: clr = FiberTracker::Normal; break;
	case 2: clr = FiberTracker::Chiral; break;
	}

	m_tracker->setFractionalAnisotropy(m_pFA);
	m_tracker->setPrincipalDirection(m_pPD);
	m_tracker->setSeedDensity(1./density);
	m_tracker->setColorFunction(clr);

	m_sceneGraph->addChild(m_tracker->track());
	m_eView->viewAll();
	m_eView->saveHomePosition();
}

void viewer3d::setSequencerSpeed(double speed) {
	m_soVolSwitch->speed = speed;
}

void viewer3d::loadPoints(const QString &filename) {
	//// read points ...
	//std::ifstream infile (filename.toAscii());
	//unsigned int temp;
	//infile.read((char *)(&temp), sizeof(unsigned int));
	//// infile.read((char *)(&temp2), sizeof(unsigned int));
	//std::cout << "Point size is " << temp << std::endl;
	//double *pts = new double[3*temp];

	//infile.read((char *)pts, sizeof(double)*3*temp);

	//infile.close();

	//// now render these ...
	//SoSeparator *pGp = new SoSeparator;

	//SoTransform *o2v = new SoTransform;
	//m_sceneGraph->addChild(pGp);
	//pGp->addChild(o2v);

	//SoCoordinate3 * coord3 = new SoCoordinate3;

	//float xyz[temp][3];

	//for (unsigned int i=0; i<temp; i++) {
	//  // std::cout << pts[3*i] << " " << pts[3*i+1] << " " << pts[3*i+2] << std::endl; 
	//  xyz[i][0] = pts[3*i];
	//  xyz[i][1] = pts[3*i+1];
	//  xyz[i][2] = pts[3*i+2];
	//}

	//coord3->point.setValues(0, temp, xyz);
	//pGp->addChild(coord3);

	//SoDrawStyle * drawstyle = new SoDrawStyle;
	//drawstyle->style.setValue(SoDrawStyle::LINES);
	//drawstyle->pointSize = 1;
	//pGp->addChild(drawstyle);

	//SoPointSet * pointset = new SoPointSet;
	//pGp->addChild(pointset);


	//// o2v->scaleFactor.setValue(0.5, 0.5, 0.5);
	//// o2v->translation.setValue(0.5, 0.5, 0.5);

	//unsigned int lng = 1<<30;
	//// o2v->translation.setValue(-0.5, -0.5, -0.5);
	//// o2v->translation.setValue(-lng/2, -lng/2, -lng/2);
	//o2v->scaleFactor.setValue(lng, lng, lng);
	////oct2Vol->translation.setValue(-lng/2, -lng/2, -lng/2);

	//pGp->addChild(new SoCube);

	// m_eView->viewAll();
	// m_eView->saveHomePosition();
}

QColor viewer3d::getBackgroundColor() {
	SbColor col = m_eView->getBackgroundColor();
	return QColor::fromRgbF(col[0], col[1], col[2]);
}


void viewer3d::setBackgroundColor(QColor col) {
	m_eView->setBackgroundColor(SbColor(col.redF(), col.greenF(), col.blueF()));
}


QList<Volume*> viewer3d::getVolumes() {
	return m_volume; 
}

void viewer3d::alignVolumes() {
	/*
	OptimizerOctPB *opt = new OptimizerOctPB();
	// OptimizerOctGD *opt = new OptimizerOctGD();
	//opt->setLearningRate(0.8);

	oct_sim *sim = new oct_sim();

	sim->setHistBits(6);
	sim->initHists();

	sim->setTargetImage(m_volume[0]);
	sim->setSourceImage(m_volume[1]);
	sim->setLowResMode(true);

	double params[6] = { 0, 0, 0, 0, 0, 0};
	double steps[6] = { 2.0, 2.0, 2.0, 0.2, 0.2, 0.2};


	opt->setSimilarityMeasure(sim);
	opt->setInitialParameters(params, steps, 6);

	opt->optimize();
	*/
	/*
	std::cout << std::endl << "Switching to high res mode" << std::endl;
	double *x = opt->getParameters ();
	opt->setInitialParameters(x, steps, 6);
	sim->setLowResMode(false);
	opt->optimize();
	*/
	/*
	double *x = opt->getParameters ();


	std::cout << "Finished Optimization" << std::endl;

	printf ("[ ");
	for (int i = 0; i < 6; i++)
	printf ("%f ", x[i]);
	printf ("]\n");


	std::cout << "Deleting sim" << std::endl;
	delete sim;
	std::cout << "Deleting opt" << std::endl;
	delete opt;
	std::cout << "All done" << std::endl;
	*/
}

void viewer3d::animateLandmarks(const QString &dirname) {
	// read in points first ...
	QString fname = dirname + "/points.txt";
	std::ifstream in(fname.toAscii());

	double fac = 0.1;

	int n, L;

	double x,y,z;

	in >> n >> L;

	float lng = (float)L;
	float sz = 1;

	// center the octree ...
	SoTransform *pts2Vol = new SoTransform;
	pts2Vol->setName("pts2Vol");

	SoSFVec3f _mn;
	SoSFVec3f _mx;

	if (m_bHaveVolume) {
		float sx, sy, sz;
		float tx, ty, tz;

		SbBox3f bbox = m_soVolumeData[0]->getVolumeSize();

		bbox.getSize( sx, sy, sz );
		sx /= lng; sy /= lng; sz /= lng;
		bbox.getOrigin( tx, ty, tz );

		_mn = bbox.getMin();
		_mx = bbox.getMax();

		pts2Vol->scaleFactor.setValue(sx, sy, sz);
		pts2Vol->translation.setValue(tx, ty, tz);
	} 

	SoSeparator *landmarks = new SoSeparator;
	SoSeparator *fixedPts = new SoSeparator;

	landmarks->addChild(pts2Vol);


	double px, py, pz;
	QList<Point> pts; 
	for (int i=0; i<n; i++) {
		// read in a point ...
		in >> x >> y >> z;

		SoSeparator * sep = new SoSeparator;

		/*
		QString str = QString::number(x) + ", " + QString::number(y) + ", " + QString::number(z);
		std::cout << "Reading point " << x << ", " << y << ", " << z << std::endl;
		QMessageBox::warning(this, "Cardiac Viewer",
		str,
		"Retry", "Quit", 0, 0, 1);
		*/
		pts.push_back(Point(x,y,z));

		float s1,s2,s3; s1=sz;
		// std::cout << "setting scale" << std::endl;
		SoTransform *trans = new SoTransform;
		//trans->scaleFactor.setValue(s1, s1, s1);
		// s1=(float)sz/2 +x; s2=(float)sz/2 + y; s3=(float)sz/2 + z;
		s1=x; s2=y; s3=z;

		// std::cout << "setting translation" << std::endl;
		trans->translation.setValue(s1, s2, s3);
		sep->addChild(trans);
		// std::cout << "adding cube under separator under level ... " << std::endl;
		SoSphere *sph = new SoSphere;
		// cube->width = sz; cube->height = sz; cube->depth = sz;
		sph->radius = sz/2;
		// std::cout << "Added cube: " << sz << " " << s1 << " " << s2 << " " << s3 << std::endl;
		sep->addChild(sph);
		fixedPts->addChild(sep);
	}
	in.close();

	landmarks->addChild(fixedPts);
	m_geoGroup->addChild(landmarks);
	// Animation 
	SoBlinker *animSwitch = new SoBlinker;
	
	SoBaseColor * col = new SoBaseColor;
	col->rgb = SbColor(0, 1, 0);

	double *fld = new double[3*(L+1)*(L+1)*(L+1)];

	/*
	QMessageBox::warning(this, "Cardiac Viewer",
		"Allocated Memory for fld",
		"Retry", "Quit", 0, 0, 1);
	*/

	char buffname[1024];
	QString fldname;
	// set up progress dialog for this ...
	QProgressDialog progress("Loading Deformation Fields...", "Abort Load", 0, 200, this);
	progress.setValue(0);
	for (unsigned int i=0; i<200; i++) {
		progress.setValue(i);
		qApp->processEvents();
		// read in the fields and deform ...
		sprintf(buffname, "Def.%.3d.raw",i);
		fldname = dirname + "/" + QString("Def.%1.raw").arg(i,3,10,QChar('0'));
		fldname = dirname + "/" + buffname;
		/*QMessageBox::warning(this, "Cardiac Viewer",
			fldname,
			"Retry", "Quit", 0, 0, 1); */
		std::ifstream fin(fldname.toAscii(), std::ios::binary);

		if (!fin.good()) {
			QMessageBox::warning(this, "Cardiac Viewer",
				fldname,
				"Retry", "Quit", 0, 0, 1);
			break;
		}

		fin.read((char *)fld, sizeof(double)*3*(L+1)*(L+1)*(L+1));
		fin.close();



		SoSeparator* frame = new SoSeparator;
		frame->addChild(col);
		for (unsigned int j=0; j<pts.size(); j++) {
			unsigned int idx = 3*( (pts[j].z()*(L+1) + pts[j].y())*(L+1) +pts[j].x() );
			SoSeparator* sep = new SoSeparator;
			px = pts[j].x() +  fac*fld[idx];
			py = pts[j].y() +  fac*fld[idx+1];
			pz = pts[j].z() +  fac*fld[idx+2];

			float s1,s2,s3;
			// std::cout << "setting scale" << std::endl;
			SoTransform *trans = new SoTransform;
			s1=px; s2=py; s3=pz;

			trans->translation.setValue(s1, s2, s3);
			sep->addChild(trans);
			SoSphere *sph = new SoSphere;
			sph->radius = sz/2;
			sep->addChild(sph);

			frame->addChild(sep);
		}
		animSwitch->addChild(frame);

	
		if (progress.wasCanceled())
			break;
	}

	progress.setValue(200);
	landmarks->addChild(animSwitch);

	m_soVolSwitch->whichChild = 0;
	animSwitch->on = TRUE;
	animSwitch->speed = 1;
	
	delete [] fld;
}