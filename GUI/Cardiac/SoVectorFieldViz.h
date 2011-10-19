//////////////////////////////////////////////////////////////////////////////
//
//  Author      : Hari Sundar
//  Date        : August 23rd, 2006
//  File        : SoVectorFieldViz.h
//  Description : Header file defining the SoVectorFieldViz class
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _SO_VECTOR_FIELD_VIZ_H_
#define _SO_VECTOR_FIELD_VIZ_H_

#include <Inventor/engines/SoSubEngine.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/fields/SoMFInt32.h>
#include <Inventor/fields/SoSFUShort.h>

#include "SFVectorField.h"
#include "SFScalarField.h"

class Point;
//////////////////////////////////////////////////////////////////////////////
//
//  Class: SoVectorFieldViz
//
//  A subclass of the SoEngine class.  This node renders a vector field either
//  using a set of gridlines deformed by the deformation field or using vector
//  glyphs that are oriented along the vector field. 
//
//////////////////////////////////////////////////////////////////////////////

class SoVectorFieldViz : public SoEngine {

  // Required Inventor macro
  SO_ENGINE_HEADER(SoVectorFieldViz);

public:
  enum RenderType { Grid2d, Grid3d, Glyph2d, Glyph3d };
	enum ColorType { MAGNITUDE, ORIENTATION, STRAIN_CIRCUMFERENTIAL, STRAIN_RADIAL };

  // field data
  SFVectorField data;
	// Original image
	SFScalarField image;

  // sparsity value to be used for rendering the 
  SoSFFloat     sparsity;
  // The render mode ... currently supported
  SoSFUShort	drawMode;
	// The coloring mode ... 
	SoSFUShort  colorMode;
  // The render Axis ... 
  SoSFUShort	gridPlane; // YZ = 0, XZ=1, XY=2
  // The render slice ...
  SoSFUShort	sliceNum;

  SoEngineOutput  points;    // (SoMFVec3f)
  // SoEngineOutput  normals;   // (SoMFVec3f) // no normals for now ...
  // colors ... per_vertex ...
  SoEngineOutput  colors;   // (SoMFInt32)
  // every entry corresponds to the number of 
  SoEngineOutput  lineSet;   // (SoMFInt32)

  // Inventor stuff
  static void   initClass          ();

                SoVectorFieldViz      ();

protected:
  virtual void  inputChanged (SoField *which);

private:
  virtual       ~SoVectorFieldViz     ();

  // Inventor method, called when the engine is to be evaluated
  virtual void  evaluate ();

  unsigned int pt2Col(Point p);
  // internal variables to hold data information instead of constantly
  // calling the getValue() routine from the class. The actual data is not
  // copied only a pointer to the data array.
  SbVec3f       min, max, voxel;
  float        *values;
  int           dims[3];
  
  // output buffer variables.  an output field cannot be set directly
  // (efficiently), it is more efficient to store everything in a separate
  // variable and write all data out at once when done.
  SoMFInt32     inds;
  SoMFInt32     cols;
  SoMFVec3f     verts;
  
};

#endif /* _SO_VECTOR_FIELD_VIZ_H_ */
