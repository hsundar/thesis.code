//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Hari Sundar
//  Date        : August 23rd, 2006
//  File        : SoVectorFieldViz.cpp
//  Description : Implementation of the SoVectorFieldViz class.
//
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstring>

#include "SoVectorFieldViz.h"
#include "Point.h"

#define BGCOLOR 0x00000000

// Inventor stuff
SO_ENGINE_SOURCE(SoVectorFieldViz);

// Inventor stuff
void
SoVectorFieldViz::initClass ()
{
  SO_ENGINE_INIT_CLASS(SoVectorFieldViz, SoEngine, "Engine");
}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Default constructor initializes all class fields
//
//////////////////////////////////////////////////////////////////////////////
SoVectorFieldViz::SoVectorFieldViz ()
{
  SO_ENGINE_CONSTRUCTOR(SoVectorFieldViz);

  dims[0] = dims[1] = dims[2] = 0;
  SO_ENGINE_ADD_INPUT(data,               (dims, 0, false));
  SO_ENGINE_ADD_INPUT(sparsity,                      (0.25));
  SO_ENGINE_ADD_INPUT(drawMode, 			(0));
	SO_ENGINE_ADD_INPUT(colorMode, 			(0));
  SO_ENGINE_ADD_INPUT(gridPlane, 			(2));
  SO_ENGINE_ADD_INPUT(sliceNum, 			(0));


  SO_ENGINE_ADD_OUTPUT(points,                   SoMFVec3f);
  // SO_ENGINE_ADD_OUTPUT(normals,                  SoMFVec3f);
  SO_ENGINE_ADD_OUTPUT(colors,                  SoMFInt32);
  SO_ENGINE_ADD_OUTPUT(lineSet,                  SoMFInt32);
  // SO_ENGINE_ADD_OUTPUT(numTriangles,             SoSFFloat);

  // The grid/glyph will be created with coordinates between 0 and 1, making it
  // easier to translate and scale
  min = SbVec3f(0.0, 0.0, 0.0);
  max = SbVec3f(1.0, 1.0, 1.0);
}

SoVectorFieldViz::~SoVectorFieldViz ()
{

}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Whenever the data field changes set the flag to recompute the
//    gradients.
//
//////////////////////////////////////////////////////////////////////////////
void
SoVectorFieldViz::inputChanged (SoField *which)
{
  // No special case here ... everything is done on the fly ...
//  if (which == &data)
//    computeGradients = true;

}

//////////////////////////////////////////////////////////////////////////////
//
//  Description:
//    Called by Inventor each time an output has been called and an input has
//    been modified.  The points and the indexes will be computed here.
//
//////////////////////////////////////////////////////////////////////////////
void
SoVectorFieldViz::evaluate ()
{

  // reset the arrays holding the vertices and indexes before they are
  // modified.
  verts.setNum(0);
  inds.setNum(0);

  // grab the current values of variables so that only one value will be used
  // during the computation.
  float spar = sparsity.getValue();
  int mode = drawMode.getValue();
	int cmode = colorMode.getValue();
  int plane = gridPlane.getValue();
  int slice = sliceNum.getValue();

  int step = (int)(1./spar);
	unsigned int col;

  values = data.getValue(dims);

	float* tissue = image.getValue(dims);

  // now based on the mode ... compute the points & lineSets
  int cnt=0, pcnt=0;

  switch (mode) {
    case 0:
      // Grid2d
      if (plane == 0) { // x axis ...
        // step through the correct slice ... in the YZ plane ...
        int i=slice;
        // first pass 
        for (int k=0; k<dims[2]; k+=step) {
          int lineCnt=0;
          for (int j=0; j<dims[1]; j+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index] + i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2] > 0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}

            cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // j
          inds.set1Value(pcnt++, lineCnt);
        } // k
        // second pass 
        for (int j=0; j<dims[1]; j+=step) {
          int lineCnt=0;
          for (int k=0; k<dims[2]; k+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index] + i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2] > 0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // k
          inds.set1Value(pcnt++, lineCnt);
        } // j

      } else if (plane == 1) { // y axis 
        // step through the correct slice ... in the XZ plane ...
        int j=slice;
        // first pass 
        for (int k=0; k<dims[2]; k+=step) {
          int lineCnt=0;
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2] > 0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // i
          inds.set1Value(pcnt++, lineCnt);
        } // k
        // second pass 
        for (int i=0; i<dims[0]; i+=step) {
          int lineCnt=0;
          for (int k=0; k<dims[2]; k+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2]>0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // k
          inds.set1Value(pcnt++, lineCnt);
        } // i
      } else { // z axis
        // step through the correct slice ... in the XY plane ...
        int k=slice;
        // first pass 
        for (int j=0; j<dims[1]; j+=step) {
          int lineCnt=0;
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2]>0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // i
          inds.set1Value(pcnt++, lineCnt);
        } // j
        // second pass 
        for (int i=0; i<dims[0]; i+=step) {
          int lineCnt=0;
          for (int j=0; j<dims[1]; j+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2]>0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          } // k
          inds.set1Value(pcnt++, lineCnt);
        } // j
      }
      break;
    case 1:
      // Grid3d -- will require 3 passes ...
      // first ... all the lines along the x axis ...
      for (int k=0; k<dims[2]; k+=step) {
        for (int j=0; j<dims[1]; j+=step) {
          int lineCnt=0;
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2]>0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          }
          inds.set1Value(pcnt++, lineCnt);
        }
      }
      // second ... all the lines along the y axis ...
      for (int k=0; k<dims[2]; k+=step) {
        for (int i=0; i<dims[0]; i+=step) {
          int lineCnt=0;
          for (int j=0; j<dims[1]; j+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2]>0.1)
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          }
          inds.set1Value(pcnt++, lineCnt);
        }
      }
      // last ... all the lines along the z axis ...
      for (int j=0; j<dims[1]; j+=step) {
        for (int i=0; i<dims[0]; i+=step) {
          int lineCnt=0;
          for (int k=0; k<dims[2]; k+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
						int index2 = i + dims[0]*(j + k*dims[1]);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
						// The color for the vertex
						if (cmode == ORIENTATION) {
							col = pt2Col(Point(values[index], values[index+1], values[index+2]));
						} else if (cmode == MAGNITUDE) {
							if (tissue[index2])
								col = 0x0000FFFF;
							else
								col = 0x00000000;
						}
						cols.set1Value(cnt, col);
            cnt++;
            lineCnt++;
          }
          inds.set1Value(pcnt++, lineCnt);
        }
      }
      
      break;
    case 2:
      // Glyph 2d
      if (plane==0) {
        int i=slice;
        for (int k=0; k<dims[2]; k+=step)
          for (int j=0; j<dims[1]; j+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
            verts.set1Value(cnt, i, j, k);
            cols.set1Value(cnt++, 0x0000FFFF);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
            cols.set1Value(cnt++, 0x0000FFFF);
            inds.set1Value(pcnt++, 2);
          }
      } else if (plane == 1) {
        int j=slice;
        for (int k=0; k<dims[2]; k+=step)
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
            verts.set1Value(cnt, i, j, k);
            cols.set1Value(cnt++, 0x0000FFFF);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
            cols.set1Value(cnt++, 0x0000FFFF);
            inds.set1Value(pcnt++, 2);
          }
      } else {
        int k=slice;
        for (int j=0; j<dims[1]; j+=step)
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
            verts.set1Value(cnt, i, j, k);
            cols.set1Value(cnt++, 0x0000FFFF);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
            cols.set1Value(cnt++, 0x0000FFFF);
            inds.set1Value(pcnt++, 2);
          }
     }
      break;
    case 3:
      // Glyph 3d
      for (int k=0; k<dims[2]; k+=step)
        for (int j=0; j<dims[1]; j+=step)
          for (int i=0; i<dims[0]; i+=step) {
            int index = 3*(i + dims[0]*(j + k*dims[1]));
            verts.set1Value(cnt, i, j, k);
            cols.set1Value(cnt++, 0x0000FFFF);
            verts.set1Value(cnt, values[index]+i, values[index+1]+j, values[index+2]+k);
            cols.set1Value(cnt++, 0x0000FFFF);
            inds.set1Value(pcnt++, 2);
          }
      break;
  }

  // send the results to the connected fields
  SO_ENGINE_OUTPUT(points, SoMFVec3f, setValues(0, verts.getNum(), verts.getValues(0)));
  // SO_ENGINE_OUTPUT(normals, SoMFVec3f, setValues(0, norms.getNum(), norms.getValues(0)));

  int i = inds.getNum();
  int j = cols.getNum();

  SO_ENGINE_OUTPUT(colors, SoMFInt32, setNum(j));
  SO_ENGINE_OUTPUT(colors, SoMFInt32, setValues(0, j, cols.getValues(0)));

  SO_ENGINE_OUTPUT(lineSet, SoMFInt32, setNum(i));
  SO_ENGINE_OUTPUT(lineSet, SoMFInt32, setValues(0, i, inds.getValues(0)));

  // SO_ENGINE_OUTPUT(numTriangles, SoSFFloat, setValue(triNum));
}

unsigned int SoVectorFieldViz::pt2Col(Point p) {
	Point cvec = p;
	cvec.normalize();
        bool trans = cvec.abs()  < 1e-3;
	cvec *= 255;

	unsigned int r, g, b, a, clr;
	r = (int)fabs(cvec.x()); g = (int)fabs(cvec.y()); b = (int)fabs(cvec.z());
        if (trans) {
          a = 0;
          clr = 0x00000000;
        }
        else {
          a = 255;
          clr = (r << 24) + (g << 16) + (b << 8) + a; 
        }
	
	// clr = (r << 24) + (g << 16) + (b << 8) + a;
	return clr;
}


